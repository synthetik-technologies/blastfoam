#!/usr/bin/env python3
"""
verify_equilibrium_chemistry.py: local chemical equilibrium in the gas-phase
chemistry (plan step C4).

Runs the single-cell closed case tests/0D/equilibriumChemistry_0D (blastFoam,
reacting nasa9/Abel-Nobel gas with 18 species, chemistry with
`equilibrium { Teq 2500; }` and no finite-rate reactions) for several initial
states, and checks, after the first time step and at the end (5 steps):

- the composition and temperature equal the equilibrium at the cell's density
  and absolute internal energy, computed with the same solver standalone
  (applications/test/Test-blastEquilibrium): |dY| < 1e-9, |dT|/T < 1e-9;
- mass, element and energy conservation: density, element abundances and
  absolute internal energy unchanged to 1e-10 (relative; energy relative to
  max(|e|, 1 MJ/kg));
- a state below Teq keeps its composition (no reactions), and equilibrium
  is applied exactly in the steps that start above Teq (the dense products
  cool below Teq in the first step).

The first-step check also covers the stage handling of the chemistry source
(the rate is evaluated at the first Runge-Kutta stage and applied in all
stages, so that the equilibrium is reached in one step).

blastEulerFoam: the 0D products-to-gas transfer case (tests/0D/
phaseTransfer_0D, compBAl products and a reacting JANAF gas in one cell) with
the gas chemistry switched on (`equilibrium { Teq 500; }`, no finite-rate
reactions). The conserved energy, sum of the phase (sensible) energies plus
the formation enthalpy of the gas minus that of the transferred products,
must stay constant to 1e-6 (relative to the sensible energy), and the gas
composition must change smoothly (the reaction rates are per unit phase
volume and are weighted by the gas volume fraction in the species equations;
without that weight the species change is 1/alpha too large and the
composition oscillates).
"""

import os
import re
import shutil
import subprocess
import sys

ROOT = "/opt/OpenFOAM/blastfoam"
HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "equilibriumChemistry_0D")
WORK = os.path.join(HERE, "equilibrium_chemistry")

TEQ = 2500

SPECIES = ["O2", "N2", "CO", "H2O", "CO2", "NH3", "CH4", "C2H6", "Al", "H2",
           "O", "N", "NO", "H", "OH", "AlO", "C(gr)", "Al2O3(s)"]

# name, T0 [K], p0 [Pa], mixture, equilibrium expected
SCENARIOS = [
    ("productsAir_4000K", 4000, 1e7, "productsAir", True),
    ("air_9000K", 9000, 1e6, "air", True),
    ("products_3000K_dense", 3000, 1e8, "products", True),
    ("air_2000K_belowTeq", 2000, 1e6, "air", False),
]


def env(cmd):
    return (
        "source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
        "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd
    )


def read_uniform(path):
    m = re.search(r"internalField\s+uniform\s+([^;]+);", open(path).read())
    return float(m.group(1))


def state(d):
    return {
        "T": read_uniform(os.path.join(d, "T")),
        "rho": read_uniform(os.path.join(d, "rho")),
        "Y": {k: read_uniform(os.path.join(d, k)) for k in SPECIES},
    }


def elements(name):
    comp = {}
    for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", name.split("(")[0]):
        comp[e] = comp.get(e, 0) + (int(n) if n else 1)
    return comp


def molweights(species_file):
    text = open(species_file).read()
    W = {}
    for m in re.finditer(r'^"?([A-Za-z0-9()]+)"?\s*\n\{(.*?)\n\}', text, re.S | re.M):
        mw = re.search(r"molWeight\s+([0-9.eE+-]+);", m.group(2))
        if mw:
            W[m.group(1)] = float(mw.group(1))
    return W


def abundances(Y, W):
    b = {}
    for k, y in Y.items():
        for e, c in elements(k).items():
            b[e] = b.get(e, 0) + c*y/W[k]
    return b


def equilibrium(case, Y, rho, T):
    """Equilibrium (and the energy of Y at T) from Test-blastEquilibrium."""
    mix = " ".join(f"({k} {v:.17g})" for k, v in Y.items() if v > 0)
    cmd = (f"Test-blastEquilibrium -file constant/thermo.nasa9Gas "
           f"-species '({' '.join(SPECIES)})' -mixture '({mix})' "
           f"-rho '({rho!r})' -Tfrozen '({T!r})'")
    out = subprocess.run(["bash", "-c", env(cmd)], capture_output=True, text=True, cwd=case)
    if out.returncode != 0:
        print(out.stdout[-2000:], out.stderr[-2000:])
        sys.exit("FAIL: Test-blastEquilibrium failed")
    lines = out.stdout.splitlines()
    i = next(j for j, l in enumerate(lines) if l.startswith("state "))
    v = lines[i].split()
    Yeq = {}
    for l in lines[i + 1:]:
        t = l.split()
        if len(t) == 2 and t[0] in SPECIES:
            Yeq[t[0]] = float(t[1])
    return {"ea": float(v[2]), "T": float(v[3]), "ok": v[5] == "true", "Y": Yeq}


def run_case(name, T0, p0, mix):
    case = os.path.join(WORK, name)
    if os.path.isdir(case):
        shutil.rmtree(case)
    shutil.copytree(CASE, case)
    r = subprocess.run(["bash", "-c", env(f"./Allrun {T0} {p0} {mix}")],
                       capture_output=True, text=True, cwd=case)
    log = open(os.path.join(case, "log.blastFoam")).read() if os.path.exists(
        os.path.join(case, "log.blastFoam")) else ""
    if r.returncode != 0 or "FOAM FATAL" in log or "\nEnd" not in log:
        print(log[-3000:])
        sys.exit(f"FAIL: {name}: blastFoam did not finish")
    times = sorted((d for d in os.listdir(case) if re.fullmatch(r"[0-9.e+-]+", d) and d != "0"),
                   key=float)
    return case, log, times


def run_euler_transfer():
    """0D blastEulerFoam products-to-gas transfer with equilibrium gas
    chemistry; returns the conserved-energy error and the largest late
    composition change between writes."""
    src = os.path.join(HERE, "phaseTransfer_0D")
    work = os.path.join(WORK, "euler")
    case = os.path.join(work, "phaseTransfer_0D")
    if os.path.isdir(work):
        shutil.rmtree(work)
    shutil.copytree(src, case, ignore=shutil.ignore_patterns("[0-9]*", "log.*", "postProcessing"))
    shutil.copy(os.path.join(HERE, "setup_case.py"), work)

    cp = os.path.join(case, "constant", "combustionProperties.gas")
    text = open(cp).read()
    open(cp, "w").write(text.replace("combustionModel  none;", "combustionModel  laminar;"))
    hdr = text.split("// * *")[0].replace("combustionProperties.gas", "chemistryProperties.gas")
    open(os.path.join(case, "constant", "chemistryProperties.gas"), "w").write(hdr + """
chemistryType
{
    solver          ode;
}

chemistry       on;

initialChemicalTimeStep 1e-08;

odeCoeffs
{
    solver          seulex;
    absTol          1e-12;
    relTol          1e-4;
}

equilibrium
{
    Teq             500;
}

reactions
{}
""")
    cd = os.path.join(case, "system", "controlDict")
    text = open(cd).read()
    open(cd, "w").write(re.sub(r"writeInterval\s+\d+;", "writeInterval   10;", text))

    r = subprocess.run(["bash", "-c", env("./Allrun")], capture_output=True, text=True, cwd=case)
    log = open(os.path.join(case, "log.blastEulerFoam")).read()
    if r.returncode != 0 or "\nEnd" not in log or "FOAM FATAL" in log:
        print(log[-3000:])
        sys.exit("FAIL: blastEulerFoam 0D transfer case did not finish")

    RU = 8314.462618
    txt = open(os.path.join(case, "constant", "thermo.compressibleGas")).read()
    species = re.search(r"species\s*\((.*?)\);", txt, re.S).group(1).split()
    hf = {}
    for sp in species:
        blk = txt[txt.index("\n" + sp + "\n"):][:2000]
        W = float(re.search(r"molWeight\s+([0-9.eE+-]+)", blk).group(1))
        a = [float(x) for x in re.search(r"lowCpCoeffs\s*\(([^)]*)\)", blk).group(1).split()]
        T = 298.15
        hf[sp] = RU/W*((((((a[4]/5*T + a[3]/4)*T + a[2]/3)*T + a[1]/2)*T + a[0])*T + a[5]))
    fv = open(os.path.join(case, "constant", "fvModels")).read()
    tf = re.search(r"targetFractions\s*\{(.*?)\}", fv, re.S).group(1)
    f = {m.group(1): float(m.group(2)) for m in re.finditer(r"(\S+)\s+([0-9.eE+-]+);", tf)}
    hprod = sum(v*hf[k] for k, v in f.items())/sum(f.values())

    times = sorted((d for d in os.listdir(case) if re.fullmatch(r"[0-9.e+-]+", d) and d != "0"),
                   key=float)
    u = lambda t, n: read_uniform(os.path.join(case, t, n))
    es = lambda t: u(t, "alphaRhoE.compBAl") + u(t, "alphaRhoE.gas")
    ech = lambda t: u(t, "alphaRho.gas")*sum(u(t, k + ".gas")*hf[k] for k in species)
    t0 = times[0]
    E0, M0 = es(t0) + ech(t0), u(t0, "alphaRho.gas")
    dE = max(abs(es(t) + ech(t) - (u(t, "alphaRho.gas") - M0)*hprod - E0)/abs(es(t0)) for t in times)
    late = times[len(times)//2:]
    dY = max(abs(u(b, k + ".gas") - u(a, k + ".gas")) for a, b in zip(late[:-1], late[1:]) for k in species)
    nEq = len(re.findall(r"Chemical equilibrium \(T > ", log))
    frozen = "frozen (not converged" in log
    return dE, dY, nEq, frozen, u(times[-1], "T.gas"), es(times[-1]) - es(t0)


def main():
    os.makedirs(WORK, exist_ok=True)
    r = subprocess.run(["bash", "-c", env("which Test-blastEquilibrium")], capture_output=True)
    if r.returncode != 0:
        subprocess.run(["bash", "-c", env(f"wmake {ROOT}/applications/test/Test-blastEquilibrium")],
                       check=True, capture_output=True)

    ok = True
    print(f"{'case':24s} {'T0':>7s} {'T':>9s} {'Teq ref':>9s} {'dT/T':>8s} {'dY':>8s} "
          f"{'mass':>8s} {'elem':>8s} {'energy':>8s}")
    for name, T0, p0, mix, expectEq in SCENARIOS:
        case, log, times = run_case(name, T0, p0, mix)
        W = molweights(os.path.join(case, "constant/thermo.nasa9Gas"))
        s0 = state(os.path.join(case, "0"))
        first = state(os.path.join(case, times[0]))
        last = state(os.path.join(case, times[-1]))
        rho = first["rho"]

        # Initial energy (the 0 directory holds the placeholder density)
        ref = equilibrium(case, s0["Y"], rho, s0["T"])
        ea0 = ref["ea"]
        scale = max(abs(ea0), 1e6)
        b0 = abundances(s0["Y"], W)

        good = True
        for label, s in (("step 1", first), ("end", last)):
            if expectEq:
                dT = abs(s["T"] - ref["T"])/ref["T"]
                dY = max(abs(s["Y"][k] - ref["Y"][k]) for k in SPECIES)
            else:
                dT = abs(s["T"] - s0["T"])/s0["T"]
                dY = max(abs(s["Y"][k] - s0["Y"][k]) for k in SPECIES)
            dm = abs(s["rho"] - rho)/rho
            b = abundances(s["Y"], W)
            de_l = max(abs(b.get(e, 0) - b0[e])/b0[e] for e in b0 if b0[e] > 0)
            ea = equilibrium(case, s["Y"], s["rho"], s["T"])["ea"]
            de = abs(ea - ea0)/scale
            g = dT < 1e-9 and dY < 1e-9 and dm < 1e-10 and de_l < 1e-10 and de < 1e-10
            if expectEq:
                g = g and ref["ok"]
            good = good and g
            print(f"{name + ' ' + label:24s} {T0:7.0f} {s['T']:9.2f} "
                  f"{(ref['T'] if expectEq else s0['T']):9.2f} {dT:8.1e} {dY:8.1e} "
                  f"{dm:8.1e} {de_l:8.1e} {de:8.1e}{'' if g else '  FAIL'}")
        # Equilibrium is applied in the steps that start above Teq
        nEq = len(re.findall(r"Chemical equilibrium \(T > ", log))
        Tstart = [s0["T"]] + [read_uniform(os.path.join(case, t, "T")) for t in times[:-1]]
        nExpected = sum(T > TEQ for T in Tstart)
        failed = "frozen (not converged" in log
        if expectEq and (nEq != nExpected or failed):
            print(f"  FAIL: equilibrium applied in {nEq} steps, expected {nExpected}"
                  + (", with failed solves" if failed else ""))
            good = False
        if not expectEq and nEq:
            print("  FAIL: equilibrium applied below Teq")
            good = False
        ok = ok and good

    dE, dY, nEq, frozen, Tgas, dEs = run_euler_transfer()
    good = dE < 1e-6 and dY < 1e-3 and nEq > 0 and not frozen
    print(f"\nblastEulerFoam 0D transfer with equilibrium gas chemistry: {nEq} equilibrium steps, "
          f"final T.gas {Tgas:.1f} K, heat release {dEs:.3e} J/m^3, conserved-energy error {dE:.1e}, "
          f"late max |dY| between writes {dY:.1e}{'' if good else '  FAIL'}")
    ok = ok and good

    print("SUCCESS: local chemical equilibrium in the gas chemistry verified" if ok
          else "FAILURE: equilibrium chemistry")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
