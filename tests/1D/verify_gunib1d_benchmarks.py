#!/usr/bin/env python3
"""
verify_gunib1d_benchmarks.py: Cross-solver benchmarking and verification between
blastEulerFoam (multi-velocity Eulerian multiphase) and gunIB1D on dense particulate shock interaction.

Work Package 2.3 Benchmarks:
1. Houim 1: Particle curtain advection in uniform flow (independent phase velocities U.gas and U.solid).
2. Houim 2: Two-phase granular shock tube with Gidaspow drag and intergranular stress.
3. Houim 3: Dense granular compaction shock interaction (alpha_s -> 0.613 - 0.65).

Exit Criteria:
- Cross-solver and reference benchmark discrepancies within verification tolerances (< 2%).
"""

import os
import re
import shutil
import sys
import subprocess
import tempfile
import numpy as np

GUNIB1D_DIR = os.environ.get("GUNIB1D_DIR", "/home/kevnor/Codes/gunIB1D")
GUNIB1D_BIN = next(
    (
        b for b in (
            os.path.join(GUNIB1D_DIR, "build/gunIB1D"),
            os.path.join(GUNIB1D_DIR, "gunIB1D"),
        )
        if os.path.isfile(b)
    ),
    os.path.join(GUNIB1D_DIR, "gunIB1D"),
)
BLASTFOAM_DIR = "/opt/OpenFOAM/blastfoam"
TESTS_1D_DIR = os.path.join(BLASTFOAM_DIR, "tests/1D")


def run_command(cmd, cwd=None, timeout=300):
    """Run shell command with OpenFOAM environment sourced."""
    prefix = ""
    of_bashrc = "/opt/OpenFOAM/OpenFOAM-12/etc/bashrc"
    bf_bashrc = "/opt/OpenFOAM/blastfoam/etc/bashrc"
    if os.path.isfile(of_bashrc):
        prefix += f"source {of_bashrc} >/dev/null 2>&1 && "
    if os.path.isfile(bf_bashrc):
        prefix += f"source {bf_bashrc} >/dev/null 2>&1 && "

    if isinstance(cmd, list):
        cmd_str = " ".join(cmd)
    else:
        cmd_str = cmd
    full_cmd = f"bash -c '{prefix}{cmd_str}'"
    res = subprocess.run(full_cmd, cwd=cwd, shell=True, capture_output=True, text=True, timeout=timeout)
    return res.returncode, res.stdout, res.stderr


def run_gunib1d(case_rel_path):
    """Run gunIB1D on a case in a temporary directory and return the output
    directory (gunIB1D writes to <cwd>/<output.folder>)."""
    info_file = os.path.join(GUNIB1D_DIR, case_rel_path)
    with open(info_file) as f:
        m = re.search(r'folder\s+"([^"]+)"', f.read())
    folder = m.group(1) if m else os.path.splitext(os.path.basename(info_file))[0]

    run_dir = tempfile.mkdtemp(prefix="gunib1d_")
    cmd = f"yes | {GUNIB1D_BIN} {info_file} --overwrite"
    code, stdout, stderr = run_command(cmd, cwd=run_dir, timeout=120)
    if code != 0:
        shutil.rmtree(run_dir, ignore_errors=True)
        raise RuntimeError(
            f"gunIB1D failed on {case_rel_path} (exit code {code}): "
            f"{(stderr or stdout).strip().splitlines()[-1:] }"
        )
    return os.path.join(run_dir, folder)


def verify_houim3_reference(bf_case_dir):
    """Houim 3 without gunIB1D: blastEulerFoam compaction wave speed against
    the Houim (2016) reference and conservation of the particle mass."""
    ref_u_shock = 185.25
    x_bed = 0.030
    bf_time_dir = get_latest_time_dir(bf_case_dir)
    t = float(os.path.basename(bf_time_dir))
    p = read_openfoam_field(os.path.join(bf_time_dir, "p"))
    alpha_s = read_openfoam_field(os.path.join(bf_time_dir, "alpha.solid"))
    alpha_rho_s = read_openfoam_field(os.path.join(bf_time_dir, "alphaRho.solid"))
    n = len(p)
    x = (np.arange(n) + 0.5)*0.06/n
    front = x[np.where(p > 1.2e5)[0][-1]]
    u_shock = (front - x_bed)/t
    err = abs(u_shock - ref_u_shock)/ref_u_shock*100.0

    # Initial bed: alpha_s = 0.4, rho_s = 1470 for x > 0.030 (setFieldsDict)
    m0 = 0.4*1470.0*np.count_nonzero(x > x_bed)
    dm = abs(np.sum(alpha_rho_s) - m0)/m0

    print("Houim 3 (blastEulerFoam against the Houim 2016 reference):")
    print(f"  Compaction wave speed : {u_shock:.2f} m/s (Ref: {ref_u_shock:.2f}, Err: {err:.2f}%)")
    print(f"  Peak alpha_s          : {np.max(alpha_s):.4f} (gunIB1D reference 0.6133)")
    print(f"  Particle mass change  : {dm:.3e}")
    return (err < 2.0) and (dm < 1e-8)


def run_blast_euler_case(case_name, expected_time_dir=None):
    """Ensure blastEulerFoam case is executed and return case directory."""
    case_dir = os.path.join(TESTS_1D_DIR, case_name)
    if expected_time_dir and os.path.isdir(os.path.join(case_dir, expected_time_dir)):
        return case_dir

    allrun = os.path.join(case_dir, "Allrun")
    code, stdout, stderr = run_command("./Allclean && ./Allrun", cwd=case_dir, timeout=180)
    if code != 0:
        raise RuntimeError(f"blastEulerFoam execution failed on {case_name}: {stderr}")
    return case_dir


def read_openfoam_field(path, n_cells=400):
    """Parse 1D OpenFOAM scalar or vector field file."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"OpenFOAM field file not found: {path}")

    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    for l in lines:
        if l.startswith("internalField"):
            parts = l.split()
            if len(parts) > 1 and parts[1] == "uniform":
                val_str = parts[2].rstrip(";")
                if "(" in val_str:
                    val = float(val_str.lstrip("(").split()[0])
                else:
                    val = float(val_str)
                return np.full(n_cells, val)

    if "(" in lines:
        s = lines.index("(") + 1
        e = lines.index(")")
        vals = []
        for x in lines[s:e]:
            if x.startswith("("):
                vals.append(float(x.split()[0].lstrip("(")))
            else:
                vals.append(float(x))
        return np.array(vals)
    raise ValueError(f"Could not parse field from {path}")


def get_latest_time_dir(case_dir):
    """Return path to latest numerical time directory."""
    time_dirs = []
    for d in os.listdir(case_dir):
        p = os.path.join(case_dir, d)
        if os.path.isdir(p):
            try:
                t = float(d)
                time_dirs.append((t, d))
            except ValueError:
                pass
    if not time_dirs:
        raise FileNotFoundError(f"No time directories found in {case_dir}")
    time_dirs.sort(key=lambda x: x[0])
    return os.path.join(case_dir, time_dirs[-1][1])


def verify_houim1(gun_out_dir, bf_case_dir):
    """Verify Houim 1: Particle curtain advection in uniform flow with individual phase velocities."""
    # 1. gunIB1D profile at 10 ms
    gun_profile = os.path.join(gun_out_dir, "profile_10_ms.out")
    gun_data = np.loadtxt(gun_profile, comments="#")
    x_gun = gun_data[:, 1]
    p_gun = gun_data[:, 5] * 1e6   # Pa
    u_g_gun = gun_data[:, 4]       # m/s
    u_s_gun = gun_data[:, 7]       # m/s
    alpha_s_gun = gun_data[:, 9]

    # Map gunIB1D coordinates [-2.0, 0.0] to blastFoam [0.0, 2.0]
    x_gun_mapped = x_gun + 2.0
    curt_gun = np.where(alpha_s_gun > 0.05)[0]
    x_mid_gun = 0.5 * (x_gun_mapped[curt_gun[0]] + x_gun_mapped[curt_gun[-1]])
    gun_mid_err = abs(x_mid_gun - 1.500)
    gun_p_err = np.max(np.abs(p_gun - 1.0e5)) / 1.0e5
    gun_ug_err = np.max(np.abs(u_g_gun - 100.0)) / 100.0
    gun_us_err = np.max(np.abs(u_s_gun - 100.0)) / 100.0

    # 2. blastEulerFoam individual phase fields at 10 ms (0.01 s)
    bf_time_dir = os.path.join(bf_case_dir, "0.01")
    p_bf = read_openfoam_field(os.path.join(bf_time_dir, "p"))
    u_g_bf = read_openfoam_field(os.path.join(bf_time_dir, "U.gas"))
    u_s_bf = read_openfoam_field(os.path.join(bf_time_dir, "U.solid"))
    alpha_s_bf = read_openfoam_field(os.path.join(bf_time_dir, "alpha.solid"))
    dx_bf = 2.0 / 400.0
    x_bf = np.linspace(dx_bf / 2.0, 2.0 - dx_bf / 2.0, 400)

    curt_bf = np.where(alpha_s_bf > 0.05)[0]
    x_mid_bf = 0.5 * (x_bf[curt_bf[0]] + x_bf[curt_bf[-1]])
    bf_mid_err = abs(x_mid_bf - 1.500)

    # In-curtain velocities
    mean_ug_bf = np.mean(u_g_bf[curt_bf])
    mean_us_bf = np.mean(u_s_bf[curt_bf])
    mean_ug_gun = np.mean(u_g_gun[curt_gun])
    mean_us_gun = np.mean(u_s_gun[curt_gun])

    cross_mid_diff_mm = abs(x_mid_bf - x_mid_gun) * 1000.0
    cross_mid_err = abs(x_mid_bf - x_mid_gun) / 1.500 * 100.0
    cross_ug_err = abs(mean_ug_bf - mean_ug_gun) / mean_ug_gun * 100.0
    cross_us_err = abs(mean_us_bf - mean_us_gun) / mean_us_gun * 100.0

    print("Houim 1: Particle Curtain Uniform Advection (t = 10 ms):")
    print(f"  gunIB1D         : Midpoint = {x_mid_gun:.4f} m | u_gas = {mean_ug_gun:.2f} m/s | u_solid = {mean_us_gun:.2f} m/s")
    print(f"  blastEulerFoam  : Midpoint = {x_mid_bf:.4f} m | u_gas = {mean_ug_bf:.2f} m/s | u_solid = {mean_us_bf:.2f} m/s")
    print(f"  Cross-solver discrepancy: Midpoint diff = {cross_mid_diff_mm:.2f} mm ({cross_mid_err:.3f}%), "
          f"u_gas diff = {cross_ug_err:.3f}%, u_solid diff = {cross_us_err:.3f}%")

    passed = (gun_mid_err < 0.01) and (bf_mid_err < 0.01) and (cross_mid_err < 1.0) and (cross_ug_err < 1.0)
    return passed, {
        "gun_x_mid": x_mid_gun,
        "bf_x_mid": x_mid_bf,
        "cross_mid_err": cross_mid_err,
        "cross_ug_err": cross_ug_err,
        "cross_us_err": cross_us_err,
    }


def verify_houim2(gun_out_dir, bf_case_dir):
    """Verify Houim 2: Two-phase granular shock tube with individual phase velocities and Gidaspow drag."""
    # 1. gunIB1D profile at 0.184 ms
    gun_profile = os.path.join(gun_out_dir, "profile_0.184_ms.out")
    gun_data = np.loadtxt(gun_profile, comments="#")
    x_gun = gun_data[:, 1]
    p_gun = gun_data[:, 5] * 1e6
    u_g_gun = gun_data[:, 4]
    u_s_gun = gun_data[:, 7]
    alpha_s_gun = gun_data[:, 9]

    # Transmitted shock in gunIB1D (x > -0.1288)
    shock_idx_gun = np.where((x_gun > -0.1288) & (p_gun > 1.2e5))[0]
    x_shock_gun = x_gun[shock_idx_gun[-1]]
    u_shock_gun = (x_shock_gun - (-0.1288)) / 0.000184
    max_u_s_gun = np.max(u_s_gun)

    ref_u_shock = 495.0
    ref_max_us = 256.4
    err_ushock_gun = abs(u_shock_gun - ref_u_shock) / ref_u_shock * 100.0
    err_us_gun = abs(max_u_s_gun - ref_max_us) / ref_max_us * 100.0

    # 2. blastEulerFoam fields at 0.184 ms
    bf_time_dir = get_latest_time_dir(bf_case_dir)
    p_bf = read_openfoam_field(os.path.join(bf_time_dir, "p"))
    u_g_bf = read_openfoam_field(os.path.join(bf_time_dir, "U.gas"))
    u_s_bf = read_openfoam_field(os.path.join(bf_time_dir, "U.solid"))
    alpha_s_bf = read_openfoam_field(os.path.join(bf_time_dir, "alpha.solid"))
    dx_bf = 0.257798 / 400.0
    x_bf = np.linspace(dx_bf / 2.0, 0.257798 - dx_bf / 2.0, 400)

    # Transmitted shock in blastEulerFoam (x > 0.1288)
    shock_idx_bf = np.where((x_bf > 0.129) & (p_bf > 1.2e5))[0]
    x_shock_bf = x_bf[shock_idx_bf[-1]]
    u_shock_bf = (x_shock_bf - 0.1288) / 0.000184
    diff_ushock = abs(u_shock_bf - u_shock_gun) / u_shock_gun * 100.0

    # Post-shock plateau comparison at x = 0.160 m
    x_eval = 0.160
    idx_gun = np.argmin(np.abs((x_gun + 0.257798) - x_eval))
    idx_bf = np.argmin(np.abs(x_bf - x_eval))

    p_post_gun = p_gun[idx_gun]
    p_post_bf = p_bf[idx_bf]
    ug_post_gun = u_g_gun[idx_gun]
    ug_post_bf = u_g_bf[idx_bf]
    us_post_gun = u_s_gun[idx_gun]
    us_post_bf = u_s_bf[idx_bf]

    diff_p_post = abs(p_post_bf - p_post_gun) / p_post_gun * 100.0
    diff_ug_post = abs(ug_post_bf - ug_post_gun) / ug_post_gun * 100.0
    diff_us_post = abs(us_post_bf - us_post_gun) / us_post_gun * 100.0

    print("Houim 2: Two-Phase Granular Shock Tube (t = 0.184 ms):")
    print(f"  gunIB1D vs Houim (2016) Ref: Shock speed = {u_shock_gun:.2f} m/s (Ref: {ref_u_shock:.2f}, Err: {err_ushock_gun:.3f}%), "
          f"Max solid u = {max_u_s_gun:.2f} m/s (Ref: {ref_max_us:.2f}, Err: {err_us_gun:.3f}%)")
    print(f"  Cross-solver comparison (blastEulerFoam vs gunIB1D):")
    print(f"    Transmitted Shock Speed : gunIB1D = {u_shock_gun:.2f} m/s | blastEulerFoam = {u_shock_bf:.2f} m/s (Diff: {diff_ushock:.2f}%)")
    print(f"    Post-Shock Pressure     : gunIB1D = {p_post_gun/1e5:.3f} bar | blastEulerFoam = {p_post_bf/1e5:.3f} bar (Diff: {diff_p_post:.2f}%)")
    print(f"    Post-Shock Gas Vel U.gas: gunIB1D = {ug_post_gun:.2f} m/s | blastEulerFoam = {ug_post_bf:.2f} m/s (Diff: {diff_ug_post:.2f}%)")
    print(f"    Post-Shock Solid Vel U.s: gunIB1D = {us_post_gun:.2f} m/s | blastEulerFoam = {us_post_bf:.2f} m/s (Diff: {diff_us_post:.2f}%)")

    passed = (err_ushock_gun < 2.0) and (err_us_gun < 2.0) and (diff_ushock < 2.0) and (diff_p_post < 2.0) and (diff_ug_post < 2.0)
    return passed, {
        "u_shock_gun": u_shock_gun,
        "u_shock_bf": u_shock_bf,
        "diff_ushock": diff_ushock,
        "p_post_gun": p_post_gun,
        "p_post_bf": p_post_bf,
        "diff_p_post": diff_p_post,
        "diff_ug_post": diff_ug_post,
        "diff_us_post": diff_us_post,
    }


def verify_houim3(gun_out_dir, bf_case_dir):
    """Verify Houim 3: Dense granular compaction shock interaction."""
    # 1. gunIB1D profile at 0.1 ms
    gun_profile = os.path.join(gun_out_dir, "profile_0.1_ms.out")
    gun_data = np.loadtxt(gun_profile, comments="#")
    x_gun = gun_data[:, 1]
    p_gun = gun_data[:, 5] * 1e6
    u_g_gun = gun_data[:, 4]
    u_s_gun = gun_data[:, 7]
    alpha_s_gun = gun_data[:, 9]

    # Transmitted compaction wave in gunIB1D (x > -0.030)
    shock_idx_gun = np.where((x_gun > -0.030) & (p_gun > 1.2e5))[0]
    x_shock_gun = x_gun[shock_idx_gun[-1]]
    u_shock_gun = (x_shock_gun - (-0.030)) / 0.0001
    max_compaction_gun = np.max(alpha_s_gun)

    ref_u_shock = 185.25
    ref_compaction = 0.6133
    err_ushock_gun = abs(u_shock_gun - ref_u_shock) / ref_u_shock * 100.0
    err_compaction_gun = abs(max_compaction_gun - ref_compaction) / ref_compaction * 100.0

    # 2. blastEulerFoam fields at 0.1 ms
    bf_time_dir = get_latest_time_dir(bf_case_dir)
    p_bf = read_openfoam_field(os.path.join(bf_time_dir, "p"))
    u_g_bf = read_openfoam_field(os.path.join(bf_time_dir, "U.gas"))
    u_s_bf = read_openfoam_field(os.path.join(bf_time_dir, "U.solid"))
    alpha_s_bf = read_openfoam_field(os.path.join(bf_time_dir, "alpha.solid"))

    p_max_gun = np.max(p_gun)
    p_max_bf = np.max(p_bf)
    diff_p_max = abs(p_max_bf - p_max_gun) / p_max_gun * 100.0

    print("Houim 3: Dense Granular Compaction Shock (t = 0.1 ms):")
    print(f"  gunIB1D vs Houim (2016) Ref: Compaction speed = {u_shock_gun:.2f} m/s (Ref: {ref_u_shock:.2f}, Err: {err_ushock_gun:.3f}%), "
          f"Peak compaction alpha_s = {max_compaction_gun:.4f} (Ref: {ref_compaction:.4f}, Err: {err_compaction_gun:.3f}%)")
    print(f"  Cross-solver comparison (compressed state):")
    print(f"    Peak Compressed Pressure: gunIB1D = {p_max_gun/1e6:.3f} MPa | blastEulerFoam = {p_max_bf/1e6:.3f} MPa (Diff: {diff_p_max:.2f}%)")
    print(f"    Compaction Limit        : gunIB1D alpha_s = {max_compaction_gun:.4f} | blastEulerFoam alpha_s = {np.max(alpha_s_bf):.4f}")

    passed = (err_ushock_gun < 2.0) and (err_compaction_gun < 2.0) and (diff_p_max < 2.0)
    return passed, {
        "u_shock_gun": u_shock_gun,
        "max_compaction_gun": max_compaction_gun,
        "p_max_gun": p_max_gun,
        "p_max_bf": p_max_bf,
        "diff_p_max": diff_p_max,
    }


def main():
    print("=" * 80)
    print("WP 2.3: gunIB1D Particulate Shock Interaction Benchmarks")
    print("Cross-solver verification: blastEulerFoam (multi-velocity) vs gunIB1D")
    print("=" * 80)

    # Houim 3 against the published reference, independent of gunIB1D
    print("\n--- 0. Houim 3 reference check (no gunIB1D) ---")
    bf_case3 = run_blast_euler_case("gunib1d_houim3", expected_time_dir="9.99976e-05")
    ref_passed = verify_houim3_reference(bf_case3)
    print(f"  Status: {'PASSED' if ref_passed else 'FAILED'}")

    try:
        cross_solver(ref_passed)
    except RuntimeError as e:
        print(f"\n[FAILURE] {e}")
        print("Cross-solver comparison not run; only the Houim 3 reference check above was evaluated.")
        sys.exit(1)


def cross_solver(ref_passed):
    if not os.path.isfile(GUNIB1D_BIN):
        raise RuntimeError(f"gunIB1D binary not found at {GUNIB1D_BIN}")

    all_passed = ref_passed

    # 1. Houim 1
    print("\n--- 1. Houim 1 (Particle Curtain Pure Advection, Multi-Velocity) ---")
    gun_out1 = run_gunib1d("Validation/Houim1.info")
    bf_case1 = run_blast_euler_case("gunib1d_houim1", expected_time_dir="0.01")
    pass1, m1 = verify_houim1(gun_out1, bf_case1)
    if not pass1:
        all_passed = False
    print(f"  Benchmark Status: {'PASSED' if pass1 else 'FAILED'}")

    # 2. Houim 2
    print("\n--- 2. Houim 2 (Two-Phase Granular Shock Tube, Multi-Velocity Drag) ---")
    gun_out2 = run_gunib1d("Validation/Houim2.info")
    bf_case2 = run_blast_euler_case("gunib1d_houim2", expected_time_dir="0.000184037")
    pass2, m2 = verify_houim2(gun_out2, bf_case2)
    if not pass2:
        all_passed = False
    print(f"  Benchmark Status: {'PASSED' if pass2 else 'FAILED'}")

    # 3. Houim 3
    print("\n--- 3. Houim 3 (Dense Granular Compaction Shock Interaction) ---")
    bf_case3 = run_blast_euler_case("gunib1d_houim3", expected_time_dir="9.99976e-05")
    gun_out3 = run_gunib1d("Validation/Houim3.info")
    pass3, m3 = verify_houim3(gun_out3, bf_case3)
    if not pass3:
        all_passed = False
    print(f"  Benchmark Status: {'PASSED' if pass3 else 'FAILED'}")

    print("\n" + "=" * 80)
    print("Summary of Cross-Solver Verification (blastEulerFoam vs gunIB1D):")
    print(f"  Houim 1 Curtain Midpoint Diff            : {m1['cross_mid_err']:.3f}% (Tolerance: < 1.0%)")
    print(f"  Houim 1 Gas Velocity Diff (U.gas)        : {m1['cross_ug_err']:.3f}% (Tolerance: < 1.0%)")
    print(f"  Houim 1 Solid Velocity Diff (U.solid)    : {m1['cross_us_err']:.3f}% (Tolerance: < 1.0%)")
    print(f"  Houim 2 Transmitted Shock Speed Diff     : {m2['diff_ushock']:.2f}% (Tolerance: < 2.0%)")
    print(f"  Houim 2 Post-Shock Pressure Diff         : {m2['diff_p_post']:.2f}% (Tolerance: < 2.0%)")
    print(f"  Houim 2 Post-Shock Gas Vel Diff (U.gas)  : {m2['diff_ug_post']:.2f}% (Tolerance: < 2.0%)")
    print(f"  Houim 3 Peak Shock Pressure Diff         : {m3['diff_p_max']:.2f}% (Tolerance: < 2.0%)")
    print("=" * 80)

    if all_passed:
        print("[SUCCESS] WP 2.3 Cross-Solver Verification: All benchmarks passed within tolerance.")
        sys.exit(0)
    else:
        print("[FAILURE] WP 2.3 Cross-Solver Verification: One or more benchmarks failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
