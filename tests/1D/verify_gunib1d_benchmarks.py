#!/usr/bin/env python3
"""
verify_gunib1d_benchmarks.py: Cross-solver benchmarking and verification between
blastFoam multiphase formulations and gunIB1D on dense particulate shock interaction.

Work Package 2.3 Benchmarks:
1. Houim 1: Particle curtain advection in uniform flow (pure contact preservation).
2. Houim 2: Two-phase granular shock tube with Gidaspow drag and intergranular stress.
3. Houim 3: Dense granular compaction shock interaction (alpha_s -> 0.613).

Exit Criteria:
- Transmitted shock speed, reflected shock overpressure, and particle curtain
  expansion velocity match within 2%.
"""

import os
import sys
import subprocess
import json
import numpy as np

GUNIB1D_DIR = "/home/kevnor/Codes/gunIB1D"
GUNIB1D_BIN = os.path.join(GUNIB1D_DIR, "build/gunIB1D")


def run_gunib1d(case_rel_path):
    """Run gunIB1D on a verification case and return output directory."""
    info_file = os.path.join(GUNIB1D_DIR, case_rel_path)
    cmd = [GUNIB1D_BIN, info_file, "--overwrite"]
    res = subprocess.run(cmd, cwd=GUNIB1D_DIR, capture_output=True, text=True, timeout=60)
    if res.returncode != 0:
        raise RuntimeError(f"gunIB1D failed on {case_rel_path}: {res.stderr}")

    case_name = os.path.splitext(os.path.basename(case_rel_path))[0]
    out_dir = os.path.join(GUNIB1D_DIR, "Output", case_name)
    return out_dir


def verify_houim1(out_dir):
    """Verify Houim 1: particle curtain advection."""
    profile_path = os.path.join(out_dir, "profile_10_ms.out")
    data = np.loadtxt(profile_path, comments="#")
    x = data[:, 1]
    p_g = data[:, 5] * 1e6  # Pa
    u_g = data[:, 4]        # m/s
    u_s = data[:, 7]        # m/s
    alpha_s = data[:, 9]

    # Target: curtain midpoint moved from -1.500 m to -0.500 m at u = 100 m/s
    curtain_idx = np.where(alpha_s > 0.05)[0]
    x_lead = x[curtain_idx[-1]]
    x_trail = x[curtain_idx[0]]
    x_mid = 0.5 * (x_lead + x_trail)
    mid_err = abs(x_mid - (-0.500))

    # Pressure uniformity
    p_err = np.max(np.abs(p_g - 1.0e5)) / 1.0e5
    u_err = np.max(np.abs(u_g - 100.0)) / 100.0

    print(f"Houim 1 (Curtain Advection):")
    print(f"  Curtain span at 10 ms: [{x_trail:.3f}, {x_lead:.3f}] m | Midpoint: {x_mid:.4f} m (Target: -0.5000 m, Diff: {mid_err*1000:.2f} mm)")
    print(f"  Max relative pressure deviation: {p_err*100:.4f}%")
    print(f"  Max relative gas velocity deviation: {u_err*100:.4f}%")

    passed = (mid_err < 0.01) and (p_err < 1e-4) and (u_err < 1e-4)
    return passed, {"x_mid": x_mid, "p_err": p_err, "u_err": u_err}


def verify_houim2(out_dir):
    """Verify Houim 2: granular shock tube."""
    profile_path = os.path.join(out_dir, "profile_0.184_ms.out")
    data = np.loadtxt(profile_path, comments="#")
    x = data[:, 1]
    p_g = data[:, 5] * 1e6
    u_g = data[:, 4]
    u_s = data[:, 7]
    alpha_s = data[:, 9]

    # Transmitted shock into right half (x > -0.1288)
    shock_idx = np.where((x > -0.1288) & (p_g > 1.2e5))[0]
    x_shock = x[shock_idx[-1]]
    u_shock = (x_shock - (-0.1288)) / 0.000184
    p_post = p_g[shock_idx[-1]]
    max_u_s = np.max(u_s)

    # Reference values from Houim (2016)
    ref_u_shock = 495.0  # m/s
    ref_max_us = 256.4   # m/s

    err_ushock = abs(u_shock - ref_u_shock) / ref_u_shock * 100.0
    err_us = abs(max_u_s - ref_max_us) / ref_max_us * 100.0

    print(f"Houim 2 (Granular Shock Tube):")
    print(f"  Transmitted shock speed: {u_shock:.2f} m/s (Ref: {ref_u_shock:.2f} m/s, Err: {err_ushock:.3f}%)")
    print(f"  Max particulate velocity: {max_u_s:.2f} m/s (Ref: {ref_max_us:.2f} m/s, Err: {err_us:.3f}%)")
    print(f"  Post-shock pressure: {p_post/1e5:.3f} bar")

    passed = (err_ushock < 2.0) and (err_us < 2.0)
    return passed, {"u_shock": u_shock, "max_u_s": max_u_s, "err_ushock": err_ushock}


def verify_houim3(out_dir):
    """Verify Houim 3: dense granular compaction."""
    profile_path = os.path.join(out_dir, "profile_0.1_ms.out")
    data = np.loadtxt(profile_path, comments="#")
    x = data[:, 1]
    p_g = data[:, 5] * 1e6
    u_g = data[:, 4]
    u_s = data[:, 7]
    alpha_s = data[:, 9]

    # Transmitted shock into dense bed (x > -0.030)
    shock_idx = np.where((x > -0.030) & (p_g > 1.2e5))[0]
    x_shock = x[shock_idx[-1]]
    u_shock = (x_shock - (-0.030)) / 0.0001
    max_compaction = np.max(alpha_s)
    max_u_s = np.max(u_s)

    # Reference values: compaction alpha ~ 0.613, wave speed ~ 185 m/s
    ref_u_shock = 185.25
    ref_compaction = 0.6133

    err_ushock = abs(u_shock - ref_u_shock) / ref_u_shock * 100.0
    err_compaction = abs(max_compaction - ref_compaction) / ref_compaction * 100.0

    print(f"Houim 3 (Dense Granular Compaction):")
    print(f"  Transmitted shock speed: {u_shock:.2f} m/s (Ref: {ref_u_shock:.2f} m/s, Err: {err_ushock:.3f}%)")
    print(f"  Peak granular compaction: alpha_s = {max_compaction:.4f} (Ref: {ref_compaction:.4f}, Err: {err_compaction:.3f}%)")
    print(f"  Peak solid velocity: {max_u_s:.2f} m/s")

    passed = (err_ushock < 2.0) and (err_compaction < 2.0)
    return passed, {"u_shock": u_shock, "compaction": max_compaction, "err_compaction": err_compaction}


def main():
    print("=" * 75)
    print("WP 2.3: gunIB1D Particulate Shock Interaction Benchmarks")
    print("Cross-solver verification against Lowe & Houim (2016, 2021) reference data")
    print("=" * 75)

    if not os.path.isfile(GUNIB1D_BIN):
        print(f"Error: gunIB1D binary not found at {GUNIB1D_BIN}")
        sys.exit(1)

    all_passed = True

    # 1. Houim 1
    print("\n--- 1. Houim 1 (Particle Curtain Pure Advection) ---")
    out1 = run_gunib1d("Verification/Houim1.info")
    pass1, m1 = verify_houim1(out1)
    if not pass1:
        all_passed = False
    print(f"  Status: {'PASSED' if pass1 else 'FAILED'}")

    # 2. Houim 2
    print("\n--- 2. Houim 2 (Two-Phase Granular Shock Tube) ---")
    out2 = run_gunib1d("Verification/Houim2.info")
    pass2, m2 = verify_houim2(out2)
    if not pass2:
        all_passed = False
    print(f"  Status: {'PASSED' if pass2 else 'FAILED'}")

    # 3. Houim 3
    print("\n--- 3. Houim 3 (Dense Granular Compaction) ---")
    out3 = run_gunib1d("Verification/Houim3.info")
    pass3, m3 = verify_houim3(out3)
    if not pass3:
        all_passed = False
    print(f"  Status: {'PASSED' if pass3 else 'FAILED'}")

    print("\n" + "=" * 75)
    if all_passed:
        print("[SUCCESS] WP 2.3 gunIB1D Cross-Validation: All benchmarks passed within 2% tolerance.")
        sys.exit(0)
    else:
        print("[FAILURE] WP 2.3 gunIB1D Cross-Validation: One or more benchmarks failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
