#!/usr/bin/env python3
"""
verify_al_flame.py: equilibrium of a stoichiometric aluminium/air mixture
(plan step D4: Al2O3 condensation in the gas chemistry).

The gas species set of the reacting cases is extended by the aluminium
sub-oxides Al2O, AlO2, Al2O2 and liquid Al2O3, Al2O3(L) (NASA data,
utilities/thermo_fitter/nasa9_species.py), next to Al, AlO and Al2O3(s).
The equilibrium solver (blastChemicalEquilibrium) treats the two condensed
oxides as pure phases and selects them by their Gibbs energy, so melting
(2327 K, 111 kJ/mol) and the condensation of the oxide vapour follow from
the data.

The mixture starts as Al(s) + air at 298.15 K (energy of Al(cr) with
Hf = 0) and is brought to equilibrium at constant volume and energy at two
densities (about 1 and 10 atm of air). The blastFoam solution
(Test-blastEquilibrium) is compared with the independent Python solution of
verify_equilibrium.py. Checks:
1. Data: Al2O3 solid/liquid Gibbs energies cross at 2327 +- 2 K; heat of
   combustion 2Al(s) + 3/2 O2 -> Al2O3(s) 31.05 +- 0.1 MJ/kg Al.
2. Equilibrium: T and p within 0.5% of the independent solution, mass
   fractions as in verify_equilibrium.py; the flame temperature is held at
   the oxide vaporisation/dissociation plateau (3500-5000 K; it rises with
   pressure) with most of
   the Al as condensed Al2O3, and is several hundred K below the
   temperature of complete combustion to Al2O3(s) without dissociation.
"""

import os
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_equilibrium as ve  # noqa: E402

WORK = os.path.join(HERE, "al_flame")
SPECIES = ["O2", "N2", "O", "N", "NO", "Al", "AlO", "Al2O", "AlO2", "Al2O2",
           "Al2O3(s)", "Al2O3(L)"]
W_AL, W_O2, W_N2 = 26.981539, 31.9988, 28.0134
Y_O2_AIR = 0.2331


def main():
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(ve.WORK, exist_ok=True)
    species_file = os.path.join(WORK, "thermo.alox")
    subprocess.run([sys.executable, os.path.join(ve.ROOT, "utilities/thermo_fitter/nasa9_species.py"),
                    "--base", ve.BASE, "--species", *SPECIES, "-o", species_file],
                   check=True, capture_output=True)
    sp = {n: ve.Species(n, d) for n, d in ve.read_species(species_file).items()}
    ok = True

    # 1. Data
    s, l = sp["Al2O3(s)"], sp["Al2O3(L)"]
    Tm = brentq(lambda T: s.g(T) - l.g(T), 1500, 3500)
    dHf = (l.h(Tm) - s.h(Tm))*ve.RU/1e6
    q = -s.h(298.15)*ve.RU/(2*W_AL)/1e6
    print(f"Al2O3 melting {Tm:.1f} K (NASA 2327 K), heat of fusion {dHf:.1f} kJ/mol; "
          f"heat of combustion {q:.3f} MJ/kg Al")
    if abs(Tm - 2327) > 2 or abs(q - 31.05) > 0.1:
        print("FAIL: Al2O3 data")
        ok = False

    # 2. Stoichiometric Al(s)/air: 4 Al + 3 O2
    m_air = 3*W_O2/Y_O2_AIR
    Y_al = 4*W_AL/(4*W_AL + m_air)
    Y_air = 1 - Y_al
    # absolute energy (blastFoam model, Al(cr) with H = 0 at 298.15 K)
    e0 = Y_air*(Y_O2_AIR*sp["O2"].e(298.15) + (1 - Y_O2_AIR)*sp["N2"].e(298.15)) \
        - Y_al*ve.RU/W_AL*298.15
    # frozen start: complete combustion to Al2O3(s) + N2 at the same energy
    Y0 = {"Al2O3(s)": Y_al + Y_air*Y_O2_AIR, "N2": Y_air*(1 - Y_O2_AIR)}
    Tc = brentq(lambda T: sum(Y*sp[n].e(T) for n, Y in Y0.items()) - e0, 300, 20000)
    print(f"Stoichiometric Al/air: Y_Al {Y_al:.4f}; complete combustion to Al2O3(s), "
          f"no dissociation: {Tc:.0f} K")

    rhos = [1.225/Y_air, 12.25/Y_air]
    bf = ve.blastfoam(species_file, SPECIES, Y0, rhos, [Tc]*len(rhos))
    refs = []
    for st in bf:
        T, p, Y = ve.python_uv(sp, SPECIES, Y0, st["vals"][0], st["vals"][1], st["vals"][2])
        refs.append({"T": T, "p": p, "Y": Y})
    good, worst, nuv = ve.compare("Al/air: blastFoam vs Python (SLSQP)", bf, refs, SPECIES)
    if not good:
        ok = False
    for st in bf:
        T_bf = st["vals"][2]
        y = st["Y"]
        cond = y["Al2O3(s)"] + y["Al2O3(L)"]
        al_cond = cond*2*W_AL/(2*W_AL + 3*15.9994)/Y_al
        print(f"  rho {st['vals'][0]:.3f} kg/m3: T {T_bf:.1f} K, p {st['vals'][3]/1e5:.2f} bar, "
              f"Al2O3(L) {y['Al2O3(L)']:.4f}, Al2O3(s) {y['Al2O3(s)']:.4f}, "
              f"Al in condensed oxide {100*al_cond:.1f}%, AlO {y['AlO']:.4f}, Al2O {y['Al2O']:.4f}, "
              f"Al {y['Al']:.4f}, O {y['O']:.4f}")
        if not (3500 <= T_bf <= 5000) or al_cond < 0.5 or Tc - T_bf < 300:
            print("FAIL: flame temperature not held by the oxide plateau")
            ok = False

    # For comparison: the species set of the cases before D4 (Al2O3(s) only)
    old = [n for n in SPECIES if n not in ("Al2O", "AlO2", "Al2O2", "Al2O3(L)")]
    for st in bf:
        T, p, Y = ve.python_uv(sp, old, Y0, st["vals"][0], st["vals"][1], st["vals"][2])
        print(f"  rho {st['vals'][0]:.3f} kg/m3 with Al2O3(s) only (no liquid, no sub-oxides): "
              f"T {T:.1f} K, p {p/1e5:.2f} bar")

    print("SUCCESS: Al/air equilibrium with condensed Al2O3" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
