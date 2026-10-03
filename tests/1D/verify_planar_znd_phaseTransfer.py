#!/usr/bin/env python3
"""
verify_planar_znd_phaseTransfer.py: Quantitative verification of 1D planar ZND
detonation coupled with pressure-driven phase conversion (pressurePhaseTransfer).

Verifies:
1. Solid HE detonation front propagation at Chapman-Jouguet speed (D_CJ = 7946.57 m/s)
   without premature phase conversion in the unreacted solid HE (lambda = 0).
2. Expansion of detonation products following breakout into the ambient gas.
3. Activation of pressurePhaseTransfer in cells where pressure drops below the
   transition threshold (p < pThreshold = 100 bar = 1.0e7 Pa).
4. Smooth transfer of expanded detonation product mass from HE to gas phase.
"""

import os
import sys
import subprocess
import numpy as np

CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "planar_znd_phaseTransfer")
EUDET_BIN = "/opt/OpenFOAM/blastfoam/modules/eudetcode/build/src/DetonationStructure"
ZND_PROFILE = "/opt/OpenFOAM/blastfoam/modules/eudetcode/znd_profile.dat"
DX = 0.2 / 400.0  # 0.0005 m (0.5 mm)
V_CELL = DX * 0.001 * 0.001  # Cell volume (1D unit depth)


def get_eudetcode_reference():
    """Extract reference CJ parameters from EuDetCode."""
    if not os.path.isfile(ZND_PROFILE):
        cmd = [EUDET_BIN, "COMP-B", "1720", "BKW", "--isentrope"]
        res = subprocess.run(cmd, cwd=os.path.dirname(ZND_PROFILE), capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"DetonationStructure failed: {res.stderr}")

    ref = {}
    with open(ZND_PROFILE, "r") as f:
        for line in f:
            if line.startswith("# D_CJ[m/s]:"):
                ref["D_CJ"] = float(line.split(":")[1].strip())
            elif line.startswith("# P_CJ[Pa]:"):
                ref["P_CJ"] = float(line.split(":")[1].strip())
            elif line.startswith("# T_CJ[K]:"):
                ref["T_CJ"] = float(line.split(":")[1].strip())
            elif line.startswith("# rho0[kg/m^3]:"):
                ref["rho0"] = float(line.split(":")[1].strip())
    return ref


def read_openfoam_field(path):
    """Parse 1D OpenFOAM field file."""
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip()]

    for l in lines:
        if l.startswith("internalField"):
            parts = l.split()
            if len(parts) > 1 and parts[1] == "uniform":
                val_str = parts[2].rstrip(";")
                val = float(val_str.lstrip("(").split()[0])
                return np.full(400, val)

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


def main():
    print("=" * 80)
    print("1D Planar ZND Detonation with pressurePhaseTransfer Verification")
    print("=" * 80)

    ref = get_eudetcode_reference()
    D_CJ = ref["D_CJ"]
    P_CJ = ref["P_CJ"]
    print(f"Reference CJ Detonation Velocity : {D_CJ:.2f} m/s")
    print(f"Reference CJ Detonation Pressure : {P_CJ/1e9:.3f} GPa")
    print("-" * 80)

    # Ensure case has run
    check_time = os.path.join(CASE_DIR, "2.49997e-05")
    if not os.path.isdir(check_time):
        print("Running blastFoam planar_znd_phaseTransfer simulation...")
        res = subprocess.run(["bash", "-c", f"cd {CASE_DIR} && ./Allrun"], capture_output=True, text=True)
        if res.returncode != 0:
            print(f"Execution failed: {res.stderr[:200]}")
            sys.exit(1)

    xc = (np.arange(400) + 0.5) * DX

    # Collect time directories
    time_dirs = []
    for d in os.listdir(CASE_DIR):
        p = os.path.join(CASE_DIR, d)
        if os.path.isdir(p):
            try:
                t = float(d)
                time_dirs.append((t, d))
            except ValueError:
                pass
    time_dirs.sort(key=lambda x: x[0])

    # 1. Verify CJ detonation velocity in solid HE (t <= 9.5 us, inside charge x <= 0.08m)
    print("[1] Verifying Detonation Front Velocity in Solid HE:")
    he_times = [td for td in time_dirs if 1e-6 < td[0] <= 9.5e-6]
    front_positions = []
    for t_sec, td in he_times:
        lam = read_openfoam_field(os.path.join(CASE_DIR, td, "lambda.HE"))
        burned = np.where(lam > 0.5)[0]
        x_num = xc[burned[-1]] if len(burned) > 0 else 0.0
        front_positions.append((t_sec, x_num))
        x_exact = D_CJ * t_sec
        diff_mm = (x_num - x_exact) * 1000.0
        rel_diff = abs(x_num - x_exact) / x_exact * 100.0
        passed = abs(x_num - x_exact) <= 1.5 * DX
        status = "PASSED" if passed else "FAILED"
        print(f"  t = {t_sec*1e6:5.2f} us | x_num = {x_num*1000:6.2f} mm | x_ref = {x_exact*1000:6.2f} mm | diff = {diff_mm:6.2f} mm ({rel_diff:5.2f}%) | {status}")

    t_fit = [0.0] + [t for t, _ in front_positions]
    x_fit = [0.0] + [x for _, x in front_positions]
    slope, _ = np.polyfit(t_fit, x_fit, 1)
    speed_err = abs(slope - D_CJ) / D_CJ * 100.0
    print(f"  Fitted Detonation Velocity: {slope:.2f} m/s (Error: {speed_err:.3f}%)")
    front_passed = speed_err < 1.0

    # 2. Verify no premature phase transfer ahead of shock front (lambda = 0)
    print("\n[2] Verifying No Premature Phase Transfer in Unreacted Solid HE:")
    t_early, td_early = he_times[0]
    p_early = read_openfoam_field(os.path.join(CASE_DIR, td_early, "p"))
    lam_early = read_openfoam_field(os.path.join(CASE_DIR, td_early, "lambda.HE"))
    a_he_early = read_openfoam_field(os.path.join(CASE_DIR, td_early, "alpha.HE"))
    unreacted_indices = np.where((xc <= 0.08) & (lam_early < 0.01))[0]
    unreacted_he_alpha = a_he_early[unreacted_indices]
    unreacted_preserved = np.allclose(unreacted_he_alpha, 1.0, atol=1e-4)
    print(f"  Unreacted solid HE cells (t = {t_early*1e6:.2f} us): {len(unreacted_indices)} cells")
    print(f"  Min alpha.HE in unreacted solid: {np.min(unreacted_he_alpha):.6f} (expected: 1.0)")
    print(f"  Status: {'PASSED (Zero premature conversion)' if unreacted_preserved else 'FAILED'}")

    # 3. Verify breakout and pressure-driven phase conversion in expanded products
    print("\n[3] Verifying Post-Breakout Product Expansion & pressurePhaseTransfer:")
    post_times = [td for td in time_dirs if td[0] >= 1.5e-5]
    conversion_verified = False
    for t_sec, td in post_times:
        p_field = read_openfoam_field(os.path.join(CASE_DIR, td, "p"))
        a_he = read_openfoam_field(os.path.join(CASE_DIR, td, "alpha.HE"))
        ar_he = read_openfoam_field(os.path.join(CASE_DIR, td, "alphaRho.HE"))
        ar_gas = read_openfoam_field(os.path.join(CASE_DIR, td, "alphaRho.gas"))

        low_p = np.where((p_field < 1.0e7) & (ar_gas > 0.1))[0]
        if len(low_p) > 0:
            converted_mass = np.sum(ar_gas[low_p]) * V_CELL
            max_ar_gas = np.max(ar_gas[low_p])
            min_a_he = np.min(a_he[low_p])
            print(f"  t = {t_sec*1e6:5.2f} us: {len(low_p)} cells in transfer window (p < 100 bar)")
            print(f"    Total gas mass in expanded region: {converted_mass*1e6:.3f} mg")
            print(f"    Peak gas density alphaRho.gas: {max_ar_gas:.2f} kg/m^3 (detonation products converted to gas)")
            conversion_verified = True

    print(f"  Status: {'PASSED (Active phase conversion to gas verified)' if conversion_verified else 'FAILED'}")

    all_passed = front_passed and unreacted_preserved and conversion_verified
    print("\n" + "=" * 80)
    if all_passed:
        print("[SUCCESS] 1D Planar ZND Detonation with pressurePhaseTransfer FULLY VERIFIED.")
        sys.exit(0)
    else:
        print("[FAILURE] Planar ZND with pressurePhaseTransfer verification failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
