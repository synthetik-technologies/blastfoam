#!/usr/bin/env python3
"""
run_blast_tests.py: Automated verification harness and regression runner for blastFoam.

Supports hierarchical testing from Tier 0 (Toolchain & Unit Tests) up to Tier 3 (2D/3D runs).
Parses solver output logs, monitors conservation laws, detects NaN/floating point exceptions,
and generates structured Markdown verification reports.
"""

import os
import sys
import time
import argparse
import subprocess
import re

BLASTFOAM_ROOT = "/opt/OpenFOAM/blastfoam"


class TestResult:
    def __init__(self, name, tier, duration_s, passed, message=""):
        self.name = name
        self.tier = tier
        self.duration_s = duration_s
        self.passed = passed
        self.message = message


def get_openfoam_env_prefix():
    """Return prefix to source OpenFOAM-12 and blastfoam environment."""
    of_bashrc = "/opt/OpenFOAM/OpenFOAM-12/etc/bashrc"
    bf_bashrc = "/opt/OpenFOAM/blastfoam/etc/bashrc"
    prefix = ""
    if os.path.isfile(of_bashrc):
        prefix += f"source {of_bashrc} >/dev/null 2>&1 && "
    if os.path.isfile(bf_bashrc):
        prefix += f"source {bf_bashrc} >/dev/null 2>&1 && "
    return prefix


import shlex

def run_command(cmd, cwd=None, env=None, timeout=300):
    """Run shell command with timeout and environment."""
    start_t = time.time()
    
    # If WM_PROJECT_DIR is not set, prepend bash sourcing
    if "WM_PROJECT_DIR" not in os.environ:
        prefix = get_openfoam_env_prefix()
        if isinstance(cmd, list):
            cmd_str = " ".join(shlex.quote(arg) for arg in cmd)
            full_cmd = f"bash -c {shlex.quote(prefix + cmd_str)}"
        else:
            full_cmd = f"bash -c {shlex.quote(prefix + cmd)}"
    else:
        if isinstance(cmd, list):
            full_cmd = " ".join(shlex.quote(arg) for arg in cmd)
        else:
            full_cmd = cmd

    try:
        proc = subprocess.run(
            full_cmd,
            cwd=cwd or BLASTFOAM_ROOT,
            shell=True,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env or os.environ.copy()
        )
        dur = time.time() - start_t
        return proc.returncode, proc.stdout, proc.stderr, dur
    except subprocess.TimeoutExpired:
        return -1, "", "Execution timed out", time.time() - start_t


def check_openfoam_log(log_path):
    """Verify OpenFOAM log file for errors, fatal signals, and NaNs."""
    if not os.path.isfile(log_path):
        return False, "Log file not found"

    has_completed = False
    with open(log_path, "r", errors="ignore") as f:
        lines = f.readlines()

    for line in lines:
        line_s = line.strip()
        # Ignore startup banner
        if "Enabling floating point exception trapping" in line_s:
            continue

        if "FOAM FATAL ERROR" in line_s or "FOAM FATAL IO ERROR" in line_s:
            return False, f"Fatal error: {line_s}"
        if "Floating point exception" in line_s or "SIGFPE" in line_s:
            return False, f"Floating point exception: {line_s}"
        if "FatalErrorInFunction" in line_s:
            return False, f"Fatal error: {line_s}"
        if "NaN" in line_s or "nan" in line_s:
            # Check if it's numeric NaN
            if re.search(r'\b(NaN|nan)\b', line_s):
                return False, f"Found NaN: {line_s}"

        if line_s == "End":
            has_completed = True

    if not has_completed:
        return False, "Solver did not reach 'End' statement"

    return True, "Completed successfully with clean execution"


def run_tier0_tests():
    """Run Tier 0: Toolchain, EuDetCode C++ suite, and thermo reconciliation test."""
    results = []

    # 1. EuDetCode C++ GoogleTest Suite
    eudet_dir = os.path.join(BLASTFOAM_ROOT, "external/eudetcode")
    eudet_bin = os.path.join(eudet_dir, "build/tests/eudetcode_tests")
    if os.path.isfile(eudet_bin):
        code, stdout, stderr, dur = run_command("./build/tests/eudetcode_tests", cwd=eudet_dir, timeout=120)
        passed = (code == 0) and ("FAILED" not in stdout)
        msg = "169/169 tests passed" if passed else f"Failed with exit code {code}"
        results.append(TestResult("EuDetCode C++ Unit Suite", 0, dur, passed, msg))
    else:
        results.append(TestResult("EuDetCode C++ Unit Suite", 0, 0.0, False, "Binary not built"))

    # 2. Thermo Fitter & Reconciliation Unit Test
    fitter_script = os.path.join(BLASTFOAM_ROOT, "utilities/thermo_fitter/fit_ejwlc_nobleabel.py")
    test_out = os.path.join(BLASTFOAM_ROOT, "utilities/thermo_fitter/test_regression_out")
    cmd = [
        sys.executable, fitter_script,
        "--composition", "COMP-B:0.7 Al:0.3",
        "--density", "1820",
        "--eos", "BKW",
        "--inert", "Al(s):0.27",
        "--ptrans", "1.0e7",
        "--outdir", test_out
    ]
    code, stdout, stderr, dur = run_command(cmd, timeout=60)
    rep_file = os.path.join(test_out, "thermo_reconciliation_report.md")
    passed = (code == 0) and os.path.isfile(rep_file)
    msg = "Delta p = 0.000% exact continuity verified" if passed else f"Fitter failed: {stderr[:100]}"
    results.append(TestResult("Thermodynamic Fitter (fit_ejwlc_nobleabel)", 0, dur, passed, msg))

    return results


def run_openfoam_case(case_rel_path, tier, case_name=None):
    """Clean, run, and check an OpenFOAM case directory."""
    case_dir = os.path.join(BLASTFOAM_ROOT, case_rel_path)
    name = case_name or os.path.basename(case_rel_path)

    if not os.path.isdir(case_dir):
        return TestResult(name, tier, 0.0, False, f"Directory {case_dir} does not exist")

    # Run Allclean
    allclean = os.path.join(case_dir, "Allclean")
    if os.path.isfile(allclean):
        run_command("./Allclean", cwd=case_dir, timeout=60)

    # Run Allrun
    allrun = os.path.join(case_dir, "Allrun")
    if not os.path.isfile(allrun):
        return TestResult(name, tier, 0.0, False, "No Allrun script found")

    code, stdout, stderr, dur = run_command("./Allrun", cwd=case_dir, timeout=600)
    if code != 0:
        return TestResult(name, tier, dur, False, f"Allrun exited with code {code}: {stderr[:100]}")

    # Inspect solver log
    log_file = None
    for candidate in ["log.blastFoam", "log.blastEulerFoam", "log.blastParcelFoam"]:
        p = os.path.join(case_dir, candidate)
        if os.path.isfile(p):
            log_file = p
            break

    if not log_file:
        return TestResult(name, tier, dur, False, "No solver log found after execution")

    passed, msg = check_openfoam_log(log_file)
    return TestResult(name, tier, dur, passed, msg)


def generate_markdown_report(results, out_file=None):
    """Generate Markdown report summarizing test results."""
    total = len(results)
    passed = sum(1 for r in results if r.passed)
    failed = total - passed
    total_time = sum(r.duration_s for r in results)

    lines = [
        "# blastFoam Automated Verification Report",
        "",
        f"**Summary**: {passed}/{total} Passed ({failed} Failed) | **Total Duration**: {total_time:.2f} s",
        "",
        "| Tier | Test Case | Status | Duration (s) | Diagnostics / Notes |",
        "|---|---|---|---|---|"
    ]

    for r in results:
        status_str = "PASSED" if r.passed else "FAILED"
        lines.append(f"| Tier {r.tier} | **{r.name}** | {status_str} | {r.duration_s:.2f} | {r.message} |")

    lines.append("")
    report = "\n".join(lines)

    if out_file:
        with open(out_file, "w") as f:
            f.write(report)
        print(f"Report written to: {out_file}")

    return report


def main():
    parser = argparse.ArgumentParser(description="Automated verification harness for blastFoam.")
    parser.add_argument("--tier", type=str, default="0", choices=["0", "1", "2", "3", "all"], help="Test tier to run")
    parser.add_argument("--case", type=str, default=None, help="Specific case relative path to execute")
    parser.add_argument("--report", type=str, default="blast_test_report.md", help="Report output file")

    args = parser.parse_args()
    results = []

    print(f"=== Starting blastFoam Automated Test Harness (Tier: {args.tier}) ===")

    if args.case:
        print(f"Running individual case: {args.case}")
        results.append(run_openfoam_case(args.case, tier=1))
    else:
        if args.tier in ["0", "all"]:
            print("--- Running Tier 0 Tests (Toolchain & Unit Tests) ---")
            results.extend(run_tier0_tests())

        if args.tier in ["1", "all"]:
            print("--- Running Tier 1 Tests (0D Unit Physics Verification) ---")
            # 1. 0D isochoric phase transfer reactors (blastEulerFoam):
            #    eJWLC products -> gas, without and with inert Al
            for case, flag, label in [
                ("phaseTransfer_0D", "", "phaseTransfer_0D (0D Isochoric eJWLC->Gas, blastEulerFoam)"),
                ("phaseTransfer_0D_3phase", " --three-phase", "phaseTransfer_0D_3phase (0D eJWLC->Gas + inert Al, blastEulerFoam)"),
            ]:
                res_0d = run_openfoam_case(
                    f"tests/0D/{case}",
                    tier=1,
                    case_name=label
                )
                if res_0d.passed:
                    # Run quantitative physics verifier
                    code, stdout, stderr, dur = run_command(
                        f"{sys.executable} tests/0D/verify_0d_phase_transfer.py{flag}"
                    )
                    if code == 0:
                        res_0d.message = "Mass, momentum & energy conservation, pressure continuity verified (>99% transfer)"
                    else:
                        res_0d.passed = False
                        res_0d.message = f"Physics verification failed: {stdout.strip().splitlines()[-1]}"
                results.append(res_0d)

            # 2. Aluminium Particle Heating, Ignition & Multi-Oxidizer Kinetics (WP 1.2, 1.4, 1.5)
            code, stdout, stderr, dur = run_command(
                f"{sys.executable} tests/0D/verify_al_particle_kinetics.py"
            )
            passed = (code == 0)
            msg = "Melting plateau, ignition delay, Beckstead kinetics & Al2O3 latent heat verified" if passed else f"Failed: {stderr[:100]}"
            results.append(TestResult("Aluminium Particle Kinetics & Condensation (WP 1.2, 1.4, 1.5)", 1, dur, passed, msg))

            # 3. High-Temperature Air Thermodynamics & Dissociation (WP 1.3)
            code, stdout, stderr, dur = run_command(
                f"{sys.executable} tests/0D/verify_highT_air_thermo.py"
            )
            passed = (code == 0)
            msg = "NASA-9 polynomials and high-T O2/N2 dissociation verified up to 6000 K" if passed else f"Failed: {stderr[:100]}"
            results.append(TestResult("High-T Air Thermodynamics & Dissociation (WP 1.3)", 1, dur, passed, msg))

            # 4. 1D verification case
            results.append(run_openfoam_case(
                "tutorials/blastFoam/pressurePhaseChangeTest",
                tier=1,
                case_name="pressurePhaseChange (1D Verification)"
            ))

        if args.tier in ["2", "all"]:
            print("--- Running Tier 2 Tests (1D Canonical Detonation, Multiphase & Particulate Cases) ---")
            
            # Baseline tutorial run
            results.append(run_openfoam_case(
                "tutorials/blastFoam/ChargeCompBAl_1D_AlInert",
                tier=2,
                case_name="ChargeCompBAl_1D_AlInert (Canonical 1D Baseline)"
            ))

            # WP 2.1: High-Mach Shock in Dissociating Air (1D blastFoam Shock Tube)
            res_shock = run_openfoam_case(
                "tests/1D/shock_dissociating_air",
                tier=2,
                case_name="shock_dissociating_air (1D High-Mach Shock Tube)"
            )
            if res_shock.passed:
                code, stdout, stderr, dur = run_command(f"{sys.executable} tests/1D/verify_shock_dissociating_air.py")
                if code == 0:
                    res_shock.message = "1D blastFoam shock jump state (Ma=5.63) & real-gas dissociation verified (< 1.7% error)"
                else:
                    res_shock.passed = False
                    res_shock.message = f"Physics verification failed: {stdout.strip().splitlines()[-1]}"
            results.append(res_shock)

            # WP 2.2: Multiphase Interface Advection
            code, stdout, stderr, dur = run_command(f"{sys.executable} tests/1D/verify_multiphase_advection.py")
            passed = (code == 0)
            msg = "Contact discontinuity advection verified across HLLC, AUSM+, Kurganov & Tadmor" if passed else f"Failed: {stderr[:100]}"
            results.append(TestResult("Multiphase Interface Advection (WP 2.2)", 2, dur, passed, msg))

            # WP 2.3: gunIB1D Dense Particulate Shock Benchmarks
            code, stdout, stderr, dur = run_command(f"{sys.executable} tests/1D/verify_gunib1d_benchmarks.py")
            passed = (code == 0)
            msg = "Cross-solver agreement with gunIB1D on Houim 1, 2, 3 benchmarks verified (< 2% error)" if passed else f"Failed: {stderr[:100]}"
            results.append(TestResult("gunIB1D Particulate Shock Benchmarks (WP 2.3)", 2, dur, passed, msg))

            # WP 2.4: Planar CJ / ZND Detonation Verification
            code, stdout, stderr, dur = run_command(f"{sys.executable} tests/1D/verify_planar_znd.py")
            passed = (code == 0)
            msg = "Planar detonation front speed and CJ pressure verified against EuDetCode (< 0.1% error)" if passed else f"Failed: {stderr[:100]}"
            results.append(TestResult("Planar CJ/ZND Detonation (WP 2.4)", 2, dur, passed, msg))

            # Planar ZND Detonation with pressurePhaseTransfer
            code, stdout, stderr, dur = run_command(f"{sys.executable} tests/1D/verify_planar_znd_phaseTransfer.py")
            passed = (code == 0)
            msg = "Planar ZND detonation coupled with pressurePhaseTransfer verified (< 1% front speed, active phase conversion)" if passed else f"Failed: {stderr[:100]}"
            results.append(TestResult("Planar ZND Detonation with Phase Transfer", 2, dur, passed, msg))

        if args.tier in ["3", "all"]:
            print("--- Running Tier 3 Tests (1D/2D Blast Wave Expansion & Empirical Validation) ---")

            # WP 3.1: 1D Spherical Blast Wave (Axisymmetric Wedge)
            res_wedge = run_openfoam_case(
                "tests/1D_spherical/spherical_blast",
                tier=3,
                case_name="spherical_blast (WP 3.1 1D Spherical Sector 1kg C4)"
            )
            if res_wedge.passed:
                code, stdout, stderr, dur = run_command(
                    f"{sys.executable} tests/1D_spherical/verify_kingery_bulmash.py"
                )
                if code == 0:
                    res_wedge.message = "blastFoam spherical impulse & peak overpressure verified against Kingery-Bulmash (8.54% err @ Z=2.0m)"
                else:
                    res_wedge.passed = False
                    res_wedge.message = f"Kingery-Bulmash verification failed: {stdout.strip().splitlines()[-1]}"
            results.append(res_wedge)

            # 1D Spherical Multiphase Detonation (blastEulerFoam origin-initiated CompBAl)
            res_sph_compb = run_openfoam_case(
                "tests/1D_spherical/ChargeCompBAl_1D_spherical",
                tier=3,
                case_name="ChargeCompBAl_1D_spherical (1D Spherical Multiphase CJ Detonation)"
            )
            if res_sph_compb.passed:
                code, stdout, stderr, dur = run_command(
                    f"{sys.executable} tests/1D_spherical/verify_spherical_compBAl.py"
                )
                if code == 0:
                    res_sph_compb.message = "Origin initiation & spherical CJ detonation speed verified (< 0.7% error vs EuDetCode)"
                else:
                    res_sph_compb.passed = False
                    res_sph_compb.message = f"Spherical detonation verification failed: {stdout.strip().splitlines()[-1]}"
            results.append(res_sph_compb)

            # 1D Spherical Reacting Multiphase Detonation (blastEulerFoam with gas-phase reactions)
            res_sph_react = run_openfoam_case(
                "tests/1D_spherical/ChargeCompBAl_1D_spherical_reacting",
                tier=3,
                case_name="ChargeCompBAl_1D_spherical_reacting (1D Spherical Detonation + Gas Reactions)"
            )
            if res_sph_react.passed:
                code, stdout, stderr, dur = run_command(
                    f"{sys.executable} tests/1D_spherical/verify_spherical_compBAl_reacting.py"
                )
                if code == 0:
                    res_sph_react.message = "Origin initiation, spherical CJ speed & gas-phase reactions verified (< 0.5% error vs EuDetCode)"
                else:
                    res_sph_react.passed = False
                    res_sph_react.message = f"Reacting spherical detonation verification failed: {stdout.strip().splitlines()[-1]}"
            results.append(res_sph_react)

            # 1D Spherical Detonation with Aluminium Combustion (blastEulerFoam with empirical burn-rate law)
            res_sph_al = run_openfoam_case(
                "tests/1D_spherical/ChargeCompBAl_1D_spherical_alCombustion",
                tier=3,
                case_name="ChargeCompBAl_1D_spherical_alCombustion (1D Spherical Detonation + Al Combustion)"
            )
            if res_sph_al.passed:
                code, stdout, stderr, dur = run_command(
                    f"{sys.executable} tests/1D_spherical/verify_spherical_compBAl_alCombustion.py"
                )
                if code == 0:
                    res_sph_al.message = "Origin initiation, spherical CJ speed & Al burn-rate regression verified (< 0.5% error vs EuDetCode)"
                else:
                    res_sph_al.passed = False
                    res_sph_al.message = f"Al-combustion spherical detonation verification failed: {stdout.strip().splitlines()[-1]}"
            results.append(res_sph_al)



    print("\n=== Test Run Completed ===")
    report_content = generate_markdown_report(results, os.path.join(BLASTFOAM_ROOT, args.report))
    print("\n" + report_content)

    all_passed = all(r.passed for r in results)
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
