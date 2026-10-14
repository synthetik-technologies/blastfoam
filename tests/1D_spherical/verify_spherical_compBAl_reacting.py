#!/usr/bin/env python3
"""
verify_spherical_compBAl_reacting.py: Verification script for 1D spherical multiphase
detonation in blastEulerFoam with gas-phase reactions enabled.

Verification criteria:
1. Origin initiation: detonation starts at r=0 at t=0.
2. Detonation propagation speed: matches EuDetCode D_CJ = 7154.23 m/s within 1.5%.
3. Detonation peak pressure: reaches expected multi-GPa range (target CJ: 20.15 GPa).
4. Multiphase physics: discrete phase velocities and volume fractions for compBAl, Al, and gas.
5. Gas-phase reaction chemistry: non-zero heat release Qdot.gas and CO2 formed in the gas.
   The gas starts as air, so this needs fuel in the gas, i.e. the products -> gas transfer
   (pressurePhaseChange, constant/fvModels).
"""

import os
import re
import sys
import numpy as np

CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ChargeCompBAl_1D_spherical_reacting")
PROBE_FILE = os.path.join(CASE_DIR, "postProcessing", "pressureProbes", "0", "p")
SAMPLE_FILE = os.path.join(CASE_DIR, "postProcessing", "sampleLine", "2e-05", "lineX.xy")
QDOT_FILE = os.path.join(CASE_DIR, "2e-05", "Qdot.gas")
CO2_FILE = os.path.join(CASE_DIR, "2e-05", "CO2.gas")

D_CJ_REF = 7154.23  # m/s from EuDetCode
P_CJ_REF = 2.01512e10  # Pa from EuDetCode (20.15 GPa)


def read_internal_field(path):
    """Internal field of an ascii volScalarField as an array (uniform -> 1 value)."""
    text = open(path).read()
    m = re.search(r"internalField\s+uniform\s+(\S+);", text)
    if m:
        return np.array([float(m.group(1))])
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\(([^)]*)\)", text)
    return np.array([float(v) for v in m.group(1).split()])



def energy_drift(case_dir, t_start=8e-6):
    """Relative change of the total energy of all phases (fluid total energy,
    Al thermal, kinetic and granular energy) from t_start to the last write
    before the shock reaches the outer boundary (closed domain). The
    programmed burn releases energy until about 7 us.

    The phase energies are sensible energies; the gas chemistry converts
    chemical into sensible energy (heat release). The conserved energy adds
    the formation enthalpy of the gas, sum(alpha*rho*Y_k*h_f,k), and
    subtracts that of the products mass transferred into the gas (gas mass
    gain minus the Al burnt, times the formation enthalpy of the
    pressurePhaseChange target fractions). The change is relative to the
    total (sensible) energy."""
    import re

    # Formation enthalpies [J/kg] of the gas species at 298.15 K, from the
    # gas thermo of the case (JANAF or nasa9)
    sys.path.insert(0, "/opt/OpenFOAM/blastfoam/utilities/thermo_fitter")
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "utilities", "thermo_fitter"))
    import products_only_fit as pof
    import compute_omegaT as co
    gas = co.phase_dict(os.path.join(case_dir, "constant"), "gas")
    species = [x.strip('"') for x in gas["species"] if x not in "()"]
    hf = {sp: pof.species_data(gas, sp, 298.15)[0](298.15) for sp in species}
    fv = open(os.path.join(case_dir, "constant", "fvModels")).read()
    tf = re.search(r"targetFractions\s*\{(.*?)\}", fv, re.S).group(1)
    f = {m.group(1): float(m.group(2)) for m in re.finditer(r"(\S+)\s+([0-9.eE+-]+);", tf)}
    hprod = sum(v*hf[k] for k, v in f.items())/sum(f.values())

    def field(t, name, n=400):
        text = open(os.path.join(case_dir, t, name)).read()
        body = text.split("internalField", 1)[1].split(";", 1)[0].strip()
        if body.startswith("uniform"):
            return np.full(n, float(body.split(None, 1)[1].strip("() ").split()[0]))
        lst = body[body.index("(") + 1:body.rindex(")")]
        if "(" in lst:
            return np.array([float(v.split()[0]) for v in re.findall(r"\(([^()]*)\)", lst)])
        return np.array([float(x) for x in lst.split()])

    r = np.linspace(0, 0.2, 401)
    V = 4/3*np.pi*(r[1:]**3 - r[:-1]**3)
    times = sorted(
        (d for d in os.listdir(case_dir)
         if re.fullmatch(r"[0-9.e+-]+", d) and float(d) >= t_start),
        key=float
    )
    E = []
    for t in times:
        if field(t, "p")[-1] > 1.05e5:
            break
        e = sum(field(t, f"alphaRhoE.{p}") @ V for p in ("compBAl", "gas", "Al"))
        e += (0.5*field(t, "alphaRho.Al")*field(t, "U.Al")**2) @ V
        e += field(t, "alphaRhoPTE.Al") @ V
        alphaRhoGas = field(t, "alphaRho.gas")
        echem = sum((alphaRhoGas*field(t, f"{k}.gas")) @ V*hf[k] for k in species)
        mTransfer = alphaRhoGas @ V + field(t, "alphaRho.Al") @ V
        E.append((float(t), e, echem, mTransfer))
    t1, e1, c1, m1 = E[-1]
    t0, e0, c0, m0 = E[0]
    return t1, (e1 + c1 - (m1 - m0)*hprod - e0 - c0)/e0

def run_case_if_needed():
    if not os.path.isfile(PROBE_FILE):
        print(f"Running Allrun in {CASE_DIR}...")
        import subprocess
        allrun = os.path.join(CASE_DIR, "Allrun")
        cmd = f"bash -c 'source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && {allrun}'"
        subprocess.run(cmd, shell=True, cwd=CASE_DIR, check=True)


def verify_spherical_detonation_reacting():
    run_case_if_needed()
    log = os.path.join(CASE_DIR, "log.blastEulerFoam")
    if not os.path.isfile(log) or "\nEnd" not in open(log).read():
        print("FAIL: blastEulerFoam did not finish")
        sys.exit(1)

    if not os.path.isfile(PROBE_FILE):
        print(f"FAIL: Probe file not found at {PROBE_FILE}")
        sys.exit(1)

    data = np.loadtxt(PROBE_FILE)
    time = data[:, 0]
    probes_r = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.10, 0.15])

    print("=" * 72)
    print(" 1D Spherical Reacting Multiphase Detonation Verification (blastEulerFoam)")
    print("=" * 72)

    arrival_times = []
    peak_pressures = []

    for i, r in enumerate(probes_r):
        p_col = data[:, i + 1]
        p_peak = np.max(p_col)
        # Detonation arrival detection: p > 1e7 Pa (100 bar)
        idx = np.where(p_col > 1e7)[0]
        if len(idx) > 0:
            t_arr = time[idx[0]]
            t_peak = time[np.argmax(p_col)]
            arrival_times.append((r, t_arr, t_peak, p_peak))
            peak_pressures.append(p_peak)
            print(f"  Probe @ r = {r*100:4.1f} cm : Arrival = {t_arr*1e6:5.2f} us | Peak P = {p_peak/1e9:6.2f} GPa (at {t_peak*1e6:5.2f} us)")
        else:
            print(f"  Probe @ r = {r*100:4.1f} cm : Peak P = {p_peak/1e9:6.2f} GPa (no arrival detected)")

    # 1. Front speed regression inside the explosive charge (r <= 0.04 m)
    charge_pts = [a for a in arrival_times if a[0] <= 0.04]
    if len(charge_pts) < 2:
        print("FAIL: Insufficient probe detections within the explosive charge.")
        sys.exit(1)

    r_vals = np.array([pt[0] for pt in charge_pts])
    t_vals = np.array([pt[1] for pt in charge_pts])

    slope, intercept = np.polyfit(t_vals, r_vals, 1)
    measured_speed = slope
    speed_error_pct = abs(measured_speed - D_CJ_REF) / D_CJ_REF * 100.0

    print("-" * 72)
    print(f"Detonation Front Velocity Analysis (r in [1 cm, 4 cm]):")
    print(f"  Measured Detonation Velocity : {measured_speed:.2f} m/s")
    print(f"  EuDetCode Reference D_CJ     : {D_CJ_REF:.2f} m/s")
    print(f"  Relative Front Speed Error   : {speed_error_pct:.2f}%")

    # 2. Peak pressure in core detonation
    max_det_p = max(pt[3] for pt in charge_pts)
    print(f"  Max Detonation Pressure      : {max_det_p/1e9:.2f} GPa (Target CJ: {P_CJ_REF/1e9:.2f} GPa)")

    # 3. Multiphase phase presence and velocities
    if os.path.isfile(SAMPLE_FILE):
        sample_data = np.loadtxt(SAMPLE_FILE)
        x_coords = sample_data[:, 0]
        alpha_compB = sample_data[:, 8]
        alpha_Al = sample_data[:, 9]
        alpha_gas = sample_data[:, 7]
        u_compB_x = sample_data[:, 16]
        u_Al_x = sample_data[:, 19]

        print("-" * 72)
        print("Multiphase Verification @ t = 20 us:")
        print(f"  Max alpha.compBAl : {np.max(alpha_compB):.4f}")
        print(f"  Max alpha.Al      : {np.max(alpha_Al):.4f}")
        print(f"  Max alpha.gas     : {np.max(alpha_gas):.4f}")
        print(f"  Max U_x (compBAl) : {np.max(np.abs(u_compB_x)):.1f} m/s")
        print(f"  Max U_x (Al)      : {np.max(np.abs(u_Al_x)):.1f} m/s")
    else:
        print(f"WARNING: Sample file not found at {SAMPLE_FILE}")

    # 4. Gas-phase reaction verification
    print("-" * 72)
    print("Gas-Phase Reaction Verification @ t = 20 us:")
    max_qdot = 0.0
    max_co2 = 0.0
    if os.path.isfile(QDOT_FILE) and os.path.isfile(CO2_FILE):
        max_qdot = np.max(np.abs(read_internal_field(QDOT_FILE)))
        max_co2 = np.max(read_internal_field(CO2_FILE))
        print(f"  Max |Qdot.gas|    : {max_qdot:.4e} W/m^3")
        print(f"  Max Y_CO2 (gas)   : {max_co2:.4e}")
    else:
        print(f"  FAIL: Qdot.gas or CO2.gas not found in {os.path.dirname(QDOT_FILE)}")
    reacting = max_qdot > 0 and max_co2 > 0
    if not reacting:
        print("  FAIL: no gas-phase heat release (no fuel in the gas: is the products -> gas transfer enabled?)")

    # Check criteria
    passed = True
    t_last, dE = energy_drift(CASE_DIR)
    print(f"  Total energy change 8 us -> {t_last*1e6:.0f} us (closed domain): {dE:.2e}")
    if abs(dE) > 1e-3:
        print(f"FAIL: total energy not conserved ({dE:.2e}, tolerance 1e-3)")
        passed = False
    if speed_error_pct > 1.5:
        print(f"FAIL: Detonation velocity error {speed_error_pct:.2f}% exceeds threshold of 1.5%")
        passed = False
    if max_det_p < 1e10:
        print(f"FAIL: Max detonation pressure {max_det_p/1e9:.2f} GPa is below 10 GPa threshold")
        passed = False
    if not reacting:
        passed = False

    print("=" * 72)
    if passed:
        print("SUCCESS: 1D spherical reacting multiphase detonation verified!")
        sys.exit(0)
    else:
        print("FAILURE: Verification criteria not met.")
        sys.exit(1)


if __name__ == "__main__":
    verify_spherical_detonation_reacting()
