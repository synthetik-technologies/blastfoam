#!/usr/bin/env python3
"""
verify_nasa9_thermo.py: NASA-9 gas thermodynamics (plan step C2).

Generates the nasa9 species file of the Comp-B/Al gas phase with
utilities/thermo_fitter/nasa9_species.py, evaluates Cp, H and S of each species
and of air with the compiled nasa9Thermo class (applications/test/
Test-nasa9Thermo, built if missing) from 200 K to 25 000 K, and compares them
with the converted coefficients (implementation) and with the source data.

Pass criteria:
- compiled class vs coefficients: within 1e-5 (output precision) for every
  species and for air;
- vs source data: within 1e-5 for species with NASA-9 data on the common
  ranges, within 3% for fitted species (condensed C(gr), Al2O3(s), and AlO).
"""

import os
import subprocess
import sys

ROOT = "/opt/OpenFOAM/blastfoam"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "utilities/thermo_fitter"))
import nasa9_species as n  # noqa: E402

SPECIES = ["O2", "N2", "CO", "H2O", "CO2", "NH3", "CH4", "C2H6", "C(gr)",
           "Al2O3(s)", "Al", "H2", "O", "N", "NO", "H", "OH", "AlO"]
FITTED = {"C(gr)", "Al2O3(s)", "AlO"}
AIR = {"N2": 0.7669, "O2": 0.2331}
TS = [200, 250, 298.15, 500, 999, 1000, 1001, 3000, 5999, 6000, 6001,
      10000, 15000, 19999, 20000, 25000]
BASE = os.path.join(HERE, "../1D_spherical/ChargeCompBAl_1D_spherical_reacting/constant/thermo.compressibleGas")
WORK = os.path.join(HERE, "nasa9_thermo")


def env(cmd):
    return (
        "source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
        "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd
    )


def isnum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def main():
    os.makedirs(WORK, exist_ok=True)
    species_file = os.path.join(WORK, "thermo.nasa9Gas")
    subprocess.run(
        [sys.executable, os.path.join(ROOT, "utilities/thermo_fitter/nasa9_species.py"),
         "--base", BASE, "--species", *SPECIES, "-o", species_file],
        check=True
    )

    app = os.path.join(ROOT, "applications/test/Test-nasa9Thermo")
    r = subprocess.run(["bash", "-c", env("which Test-nasa9Thermo")], capture_output=True, text=True)
    if r.returncode != 0:
        subprocess.run(["bash", "-c", env(f"wmake {app}")], check=True, capture_output=True)

    cmd = (f"Test-nasa9Thermo -file {species_file} -species '({' '.join(SPECIES)})' "
           f"-mixture '({' '.join(f'({k} {v})' for k, v in AIR.items())})' "
           f"-T '({' '.join(str(t) for t in TS)})'")
    out = subprocess.run(["bash", "-c", env(cmd)], capture_output=True, text=True, cwd=WORK)
    if out.returncode != 0:
        print(out.stderr[-2000:])
        print("FAIL: Test-nasa9Thermo failed")
        sys.exit(1)
    rows = [r for r in (l.split() for l in out.stdout.splitlines())
            if len(r) == 5 and all(isnum(x) for x in r[1:])]

    eu = n.load_eudet(os.path.join(ROOT, "external/eudetcode/Data/GasFlame.json"))
    ca = n.load_cantera(["/usr/local/share/cantera/data/nasa_gas.yaml",
                         "/usr/local/share/cantera/data/nasa_condensed.yaml"])
    conv = {}
    for name in SPECIES:
        src, origin = n.source(name, eu, ca)
        c, _ = n.on_common_ranges(src)
        conv[name] = (src, n.Piecewise([(lo, hi, a) for lo, hi, a in zip(n.T_RANGES[:-1], n.T_RANGES[1:], c)]), origin)

    def rel(got, ref, T):
        return max(abs(got[0] - ref.cp(T))/abs(ref.cp(T)),
                   abs(got[1] - ref.h(T))/max(abs(ref.h(T)), 1e3),
                   abs(got[2] - ref.s(T))/max(abs(ref.s(T)), 1e-3))

    worst_c, worst_s = {}, {}
    for name, T, cp, h, s in rows:
        if name == "mixture":
            continue
        T = float(T)
        got = (float(cp), float(h), float(s))
        src, fit, _ = conv[name]
        worst_c[name] = max(worst_c.get(name, 0), rel(got, fit, T))
        worst_s[name] = max(worst_s.get(name, 0), rel(got, src, T))

    Ru = 8314.462618
    R = {k: Ru/n.molweight(k) for k in AIR}
    Rm = sum(AIR[k]*R[k] for k in AIR)
    worst_air = 0
    for name, T, cp, h, s in rows:
        if name != "mixture":
            continue
        T = float(T)
        ref = [sum(AIR[k]*R[k]*getattr(conv[k][1], q)(T) for k in AIR)/Rm for q in ("cp", "h", "s")]
        got = [float(cp), float(h), float(s)]
        worst_air = max(worst_air, max(abs(g - r)/max(abs(r), 1e3 if i == 1 else 1e-30)
                                       for i, (g, r) in enumerate(zip(got, ref))))

    ok = True
    print(f"{'species':9s} {'vs coeffs':>10s} {'vs source':>10s}  source")
    for name in SPECIES:
        tol = 3e-2 if name in FITTED else 1e-5
        good = worst_c[name] < 1e-5 and worst_s[name] < tol
        ok = ok and good
        print(f"{name:9s} {worst_c[name]:10.1e} {worst_s[name]:10.1e}  {conv[name][2]}{'' if good else '  FAIL'}")
    print(f"{'air':9s} {worst_air:10.1e}")
    ok = ok and worst_air < 1e-5
    print("SUCCESS: nasa9 thermodynamics verified to 25 000 K" if ok else "FAILURE: nasa9 thermodynamics")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
