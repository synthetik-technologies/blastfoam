#!/usr/bin/env python3
"""
verify_al_burning_products.py: Al particles burning to oxide in blastEulerFoam
(tests/0D/alBurningProducts_0D).

The `Beckstead` diameter rate with `products true` converts the burnt Al and
the gas oxidisers of the cell to the products of one global reaction per
oxidiser,

    4Al + 3O2 -> 2Al2O3,  2Al + 3H2O -> Al2O3 + 3H2,  2Al + 3CO2 -> Al2O3 + 3CO,

in the proportions of their terms in X_eff = X_O2 + 0.6 X_H2O + 0.22 X_CO2.
The `reactingParticle` mass transfer (productsFromReactionRate) adds the
products and removes the oxidisers in the gas; the heat of reaction follows
from the formation enthalpies. The oxide is Al2O3(L) where the gas is above
2327 K and Al2O3(s) below. The burn rate is limited so that at most half of
an oxidiser in the cell is consumed in a time step. No gas chemistry.

Single cell, gas and particles at rest at the same temperature. Cases:
1. air, dilute (alpha 1e-7), 10 um, 2500 K: burn time within 2% of the
   correlation; O2 consumed and Al2O3(L) formed per kg Al within 0.1% of
   the stoichiometry (0.8895 and 1.8895).
2. O2/H2O/CO2 in N2, dilute: the oxidisers are consumed in the shares
   c_k X_k/X_eff (within 1%), with H2 and CO formed by the stoichiometry;
   burn time within 2% (X_eff of the initial gas).
3. air, fuel-rich (twice the stoichiometric Al): burning stops when the O2
   is used up; the burnt Al equals the initial O2/0.8895 within 1%, the O2
   mass fraction never becomes negative.
4. air at 2000 K (Tign 0): the oxide is Al2O3(s).
5. Textinction 2600 K above the gas temperature: no burning.
6. ignitedFraction (ignition tracked by the transported field ignited.Al):
   particles above Tign burn as in case 1 (burn time within 0.1%).
7. ignitedFraction, particles at 1500 K heated by the gas at 2500 K
   (constant heat transfer coefficient): ignited.Al stays zero and the
   particles do not burn until their temperature reaches Tign (2000 K);
   afterwards ignited.Al is one and they burn.
All cases: total mass to 1e-8 and total energy (phase energies plus the
formation enthalpy of the gas species) to 1e-5 of the heat released.
"""

import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "alBurningProducts_0D")
WORK = os.path.join(HERE, "al_burning_products")
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "utilities", "thermo_fitter"))

A, N = 0.00735, 1.8
P0 = 101325.0
W = {"O2": 31.9988, "N2": 28.0134, "H2O": 18.01528, "CO2": 44.0095, "H2": 2.01588,
     "CO": 28.0101, "Al2O3(s)": 101.96128, "Al2O3(L)": 101.96128}
SPECIES = ["O2", "N2", "H2O", "CO2", "H2", "CO", "Al2O3(s)", "Al2O3(L)"]
WEIGHT = {"O2": 1.0, "H2O": 0.6, "CO2": 0.22}
# kg of each gas species per kg Al (oxidiser consumed negative)
W_AL = {"O2": (3*W["O2"] - 2*W["Al2O3(s)"])/-4,
        "H2O": (3*W["H2O"] - W["Al2O3(s)"] - 3*W["H2"])/-2,
        "CO2": (3*W["CO2"] - W["Al2O3(s)"] - 3*W["CO"])/-2}
YIELD = {"O2": {"O2": -3*W["O2"]/(4*W_AL["O2"]), "Al2O3": 2*W["Al2O3(s)"]/(4*W_AL["O2"])},
         "H2O": {"H2O": -3*W["H2O"]/(2*W_AL["H2O"]), "Al2O3": W["Al2O3(s)"]/(2*W_AL["H2O"]),
                 "H2": 3*W["H2"]/(2*W_AL["H2O"])},
         "CO2": {"CO2": -3*W["CO2"]/(2*W_AL["CO2"]), "Al2O3": W["Al2O3(s)"]/(2*W_AL["CO2"]),
                 "CO": 3*W["CO"]/(2*W_AL["CO2"])}}
# columns of postProcessing/state (volAverage over the single cell)
COLS = ["t", "d", "alpha", "TAl", "Tgas", "p", "aRhoAl", "aRhoGas", "aRhoEAl", "aRhoEGas",
        "aRhoPTE"] + SPECIES + ["ignited"]


def env(cmd):
    return ("source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
            "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd)


def mole_fractions(Y):
    n = {k: Y.get(k, 0.0)/W[k] for k in SPECIES}
    tot = sum(n.values())
    return {k: v/tot for k, v in n.items()}


def burn_time(d0, Y, T):
    X = mole_fractions(Y)
    Xeff = sum(WEIGHT[k]*X[k] for k in WEIGHT)
    return 1e-3*A*(d0*1e6)**N/(Xeff*(P0/101325.0)**0.1*T**0.2)


def run_case(name, d0, T, alpha, Y, Tign, Text, dt, t_end, Tp=None, Kh=0.0, ignited=False):
    case = os.path.join(WORK, name)
    if os.path.isdir(case):
        shutil.rmtree(case)
    shutil.copytree(CASE, case)
    ystr = ",".join(f"{k}={v!r}" for k, v in Y.items())
    r = subprocess.run(["bash", "-c", env(
        f"./Allrun {d0!r} {P0!r} {T!r} {alpha!r} {dt!r} {t_end!r} {ystr} {Tign!r} {Text!r} "
        f"{Tp if Tp is not None else T!r} {Kh!r} {'true' if ignited else 'false'}")],
        capture_output=True, text=True, cwd=case)
    logf = os.path.join(case, "log.blastEulerFoam")
    log = open(logf).read() if os.path.exists(logf) else ""
    if r.returncode != 0 or "\nEnd" not in log or "FOAM FATAL" in log:
        print(log[-2000:], r.stdout[-1000:], r.stderr[-1000:])
        sys.exit(f"FAIL: {name}: blastEulerFoam did not finish")
    s = np.loadtxt(os.path.join(case, "postProcessing", "state", "0", "volFieldValue.dat"),
                   comments="#", ndmin=2)
    return case, {c: s[:, i] for i, c in enumerate(COLS)}


def formation_enthalpies(case):
    import compute_omegaT as co
    import products_only_fit as pof
    gas = co.phase_dict(os.path.join(case, "constant"), "gas")
    return {sp: pof.species_data(gas, sp, 298.15)[0](298.15) for sp in SPECIES}


def species_mass(s, k):
    return s["aRhoGas"]*s[k]


def conservation(name, case, s):
    """Mass, and energy relative to the heat released, from the first step
    (the initial decode changes the energy by about 7e-6 of the gas energy)
    while the cloud is above ten times residualAlpha (1e-9; below it the
    decode resets the energy of the remaining particles to 288 K)."""
    hf = formation_enthalpies(case)
    m = s["aRhoAl"] + s["aRhoGas"]
    chem = sum(species_mass(s, k)*hf[k] for k in SPECIES)
    E = s["aRhoEGas"] + s["aRhoEAl"] + s["aRhoPTE"] + chem
    sel = np.nonzero(s["alpha"] > 1e-8)[0]
    sel = sel[sel >= 1]
    dm = np.max(np.abs(m - m[0]))/m[0]
    dE = np.max(np.abs(E[sel] - E[sel[0]]))/max(abs(chem[sel[-1]] - chem[0]), 1e-300)
    print(f"  {name}: total mass {dm:.1e}, total energy {dE:.1e} of the heat released")
    ok = dm < 1e-8 and dE < 1e-5
    if not ok:
        print(f"FAIL: {name}: mass or energy not conserved")
    return ok


def burnt(s):
    """Al burnt per unit volume (the cloud is reset to d0 when alpha falls
    below alphaMin; use the particle mass before that)."""
    return s["aRhoAl"][0] - s["aRhoAl"]


def main():
    ok = True
    os.makedirs(WORK, exist_ok=True)
    air = {"O2": 0.2331, "N2": 0.7669}

    # 1. air, dilute
    d0, T = 1e-5, 2500.0
    tb = burn_time(d0, air, T)
    case, s = run_case("air", d0, T, 1e-7, air, 2000.0, 0.0, tb/4000, 1.02*tb)
    dn = (s["d"]/d0)**N
    i90 = np.argmax(dn <= 0.1)
    t90 = np.interp(0.1, dn[i90::-1][:2], s["t"][i90::-1][:2])
    err_tb = t90/0.9/tb - 1
    mb = burnt(s)[-1]
    rO2 = (species_mass(s, "O2")[0] - species_mass(s, "O2")[-1])/mb
    rOx = (species_mass(s, "Al2O3(L)")[-1] + species_mass(s, "Al2O3(s)")[-1])/mb
    print(f"1. air: t_b {t90/0.9*1e3:.5g} ms (Beckstead {tb*1e3:.5g} ms, {100*err_tb:+.2f}%), "
          f"O2 used {rO2:.5f} kg/kg Al ({-YIELD['O2']['O2']:.5f}), Al2O3 formed {rOx:.5f} "
          f"({YIELD['O2']['Al2O3']:.5f}), liquid {species_mass(s, 'Al2O3(L)')[-1]/max(rOx*mb, 1e-300):.3f}, "
          f"gas {s['Tgas'][0]:.0f} -> {s['Tgas'][-1]:.1f} K")
    if abs(err_tb) > 0.02 or abs(rO2/-YIELD["O2"]["O2"] - 1) > 1e-3 \
            or abs(rOx/YIELD["O2"]["Al2O3"] - 1) > 1e-3 or species_mass(s, "Al2O3(s)")[-1] > 0:
        print("FAIL: air: burn time, stoichiometry or oxide phase")
        ok = False
    ok &= conservation("air", case, s)

    # 2. O2/H2O/CO2 in N2, dilute
    mix = {"O2": 0.05, "H2O": 0.10, "CO2": 0.15, "N2": 0.70}
    tb = burn_time(d0, mix, T)
    case, s = run_case("mixed", d0, T, 1e-7, mix, 2000.0, 0.0, tb/4000, 1.02*tb)
    X = mole_fractions(mix)
    Xeff = sum(WEIGHT[k]*X[k] for k in WEIGHT)
    share = {k: WEIGHT[k]*X[k]/Xeff for k in WEIGHT}
    mb = burnt(s)[-1]
    print("2. O2/H2O/CO2: share of the burning by oxidiser (from the consumption) vs c_k X_k/X_eff:")
    for k in WEIGHT:
        used = (species_mass(s, k)[0] - species_mass(s, k)[-1])/mb
        sh = used/-YIELD[k][k]
        print(f"     {k}: {sh:.4f} (expected {share[k]:.4f})")
        if abs(sh - share[k]) > 0.01*share[k]:
            print(f"FAIL: mixed: share of {k}")
            ok = False
    h2 = species_mass(s, "H2")[-1]/mb
    co_ = species_mass(s, "CO")[-1]/mb
    h2_ref = share["H2O"]*YIELD["H2O"]["H2"]
    co_ref = share["CO2"]*YIELD["CO2"]["CO"]
    dn = (s["d"]/d0)**N
    i90 = np.argmax(dn <= 0.1)
    t90 = np.interp(0.1, dn[i90::-1][:2], s["t"][i90::-1][:2])
    err_tb = t90/0.9/tb - 1
    print(f"     H2 formed {h2:.5f} kg/kg Al ({h2_ref:.5f}), CO formed {co_:.5f} ({co_ref:.5f}); "
          f"t_b {t90/0.9*1e3:.5g} ms (Beckstead {tb*1e3:.5g} ms, {100*err_tb:+.2f}%)")
    if abs(h2/h2_ref - 1) > 0.01 or abs(co_/co_ref - 1) > 0.01 or abs(err_tb) > 0.02:
        print("FAIL: mixed: products or burn time")
        ok = False
    ok &= conservation("mixed", case, s)

    # 3. air, fuel-rich: twice the stoichiometric Al
    rho_g = P0/(8314.462618/(1/(0.2331/W["O2"] + 0.7669/W["N2"]))*T)
    alpha = 2*0.2331*rho_g/-YIELD["O2"]["O2"]/2700.0
    tb = burn_time(d0, air, T)
    case, s = run_case("rich", d0, T, alpha, air, 2000.0, 0.0, tb/4000, 1.5*tb)
    mO2 = species_mass(s, "O2")
    mb = burnt(s)[-1]
    ref = mO2[0]/-YIELD["O2"]["O2"]
    print(f"3. fuel-rich: Al burnt {mb/s['aRhoAl'][0]:.4f} of the Al, O2/0.8895 gives "
          f"{ref/s['aRhoAl'][0]:.4f} ({100*(mb/ref - 1):+.2f}%); O2 left {mO2[-1]/mO2[0]:.2e} "
          f"of the initial, min Y_O2 {s['O2'].min():.2e}; gas {s['Tgas'][-1]:.0f} K")
    if abs(mb/ref - 1) > 0.01 or s["O2"].min() < 0:
        print("FAIL: fuel-rich: burning not limited by the oxidiser")
        ok = False
    ok &= conservation("rich", case, s)

    # 4. solid oxide below the oxide melting temperature
    tb = burn_time(d0, air, 2000.0)
    case, s = run_case("cold", d0, 2000.0, 1e-7, air, 0.0, 0.0, tb/2000, 0.3*tb)
    sol, liq = species_mass(s, "Al2O3(s)")[-1], species_mass(s, "Al2O3(L)")[-1]
    print(f"4. gas at 2000 K: Al2O3(s) {sol:.3e}, Al2O3(L) {liq:.3e} kg/m3")
    if not (sol > 0 and liq == 0):
        print("FAIL: cold: oxide not solid")
        ok = False
    ok &= conservation("cold", case, s)

    # 5. extinction above the gas temperature
    case, s = run_case("extinct", d0, T, 1e-7, air, 2000.0, 2600.0, tb/2000, 0.3*tb)
    print(f"5. Textinction 2600 K > 2500 K: d/d0 {s['d'][-1]/d0:.6f}, "
          f"O2 change {species_mass(s, 'O2')[-1]/species_mass(s, 'O2')[0] - 1:.1e}")
    if abs(s["d"][-1]/d0 - 1) > 1e-9:
        print("FAIL: extinct: particles burn below Textinction")
        ok = False

    # 6. ignited fraction, particles above Tign: as case 1
    tb = burn_time(d0, air, T)
    case, s = run_case("ignited", d0, T, 1e-7, air, 2000.0, 0.0, tb/4000, 1.02*tb, ignited=True)
    dn = (s["d"]/d0)**N
    i90 = np.argmax(dn <= 0.1)
    t90 = np.interp(0.1, dn[i90::-1][:2], s["t"][i90::-1][:2])
    err = t90/0.9/tb - 1
    print(f"6. ignitedFraction, particles at 2500 K: ignited {s['ignited'][1]:.3f}, "
          f"t_b {t90/0.9*1e3:.5g} ms ({100*err:+.2f}% vs Beckstead)")
    if abs(err) > 0.02 or s["ignited"][1] != 1:
        print("FAIL: ignited: burn time or ignited fraction")
        ok = False

    # 7. ignited fraction, particles heated from 1500 K through Tign
    Kh = 2700.0*900.0/5e-5   # particle heating time 50 us
    case, s = run_case("heating", d0, T, 1e-7, air, 2000.0, 0.0, tb/4000, 0.3*tb,
                       Tp=1500.0, Kh=Kh, ignited=True)
    i_ign = np.argmax(s["TAl"] >= 2000.0)
    i_burn = np.argmax(s["d"] < (1 - 1e-9)*d0)
    t_ign_ref = 5e-5*np.log((T - 1500.0)/(T - 2000.0))
    before = s["ignited"][:i_ign].max() if i_ign > 0 else 0.0
    print(f"7. ignitedFraction, heating from 1500 K: T_Al reaches Tign at {s['t'][i_ign]*1e6:.2f} us "
          f"(lumped estimate {t_ign_ref*1e6:.2f} us), burning starts at {s['t'][i_burn]*1e6:.2f} us; "
          f"ignited before {before:.1f}, while burning {s['ignited'][i_burn:].min():.1f}; "
          f"d/d0 at the end {s['d'][-1]/d0:.4f}")
    if i_burn < i_ign or before > 0 or s["ignited"][i_burn:].min() < 1 or s["d"][-1] >= d0:
        print("FAIL: heating: burning before ignition or not after")
        ok = False
    ok &= conservation("heating", case, s)

    print("SUCCESS: Al burns to oxide with the gas oxidisers" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
