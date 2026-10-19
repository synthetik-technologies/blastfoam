#!/usr/bin/env python3
"""
setup_e2.py: set up a copy of this case for plan step E2 (inert vs reacting
Comp-B/Al with the aluminium models of Phase D).

    python3 setup_e2.py --variant inert|products|full
        [--pThreshold 1e7] [--tau 2e-5] [--Tign 2050]

Run on a copy of the case (it edits constant/ in place). All variants use:
- Al: eMelting thermodynamics (solid/liquid Cv, latent heat at 933.47 K),
  Loth drag and RanzMarshall heat transfer with the Kavanau Nusselt number
  for (Al in gas), Beckstead burning (Tign 2050 K, burning limited to the gas
  fraction of the surrounding fluid) to Al2O3 with the gas oxidisers O2, H2O
  and CO2 (products true; the mass transfer adds Al2O3, H2 and CO and removes
  the oxidisers, with the heat of reaction), ignition tracked by the
  transported ignited fraction (ignitedFraction), with alphaMin 1e-6, and an
  alphaCutoff of 1e-6 instead of 1e-2 for the Al/gas drag (with 1e-2 the
  dilute cloud has no drag);
- gas: the 24 species of the case plus the Al sub-oxides Al2O, AlO2, Al2O2
  and liquid Al2O3 (nasa9, equilibrium above Teq).
Variants:
- inert:    no products -> gas transfer, Al does not ignite (Tign 1e5 K);
- products: transfer, Al does not ignite;
- full:     transfer and Al burning.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

CASE = os.getcwd()
BLAST = os.environ.get("BLAST_DIR", "/opt/OpenFOAM/blastfoam")
EXTRA_SPECIES = ["Al2O", "AlO2", "Al2O2", "Al2O3(L)"]


def edit(path, old, new, count=1):
    with open(path) as f:
        s = f.read()
    if s.count(old) != count:
        sys.exit(f"{path}: expected {count} occurrence(s) of:\n{old}")
    with open(path, "w") as f:
        f.write(s.replace(old, new))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--variant", choices=("inert", "products", "full"), required=True)
    ap.add_argument("--pThreshold", type=float, default=None)
    ap.add_argument("--tau", type=float, default=None)
    ap.add_argument("--Tign", type=float, default=2050.0,
                    help="ignition temperature of the full variant [K]")
    a = ap.parse_args()
    c = os.path.join(CASE, "constant")

    # Al thermodynamics
    edit(os.path.join(c, "physicalProperties.Al"), "    thermo          eConst;", "    thermo          eMelting;")
    edit(os.path.join(c, "physicalProperties.Al"),
         "    type            eConst;\n    Cv              900;\n    Hf              0;",
         "    // NASA (CODATA) Al(cr)/Al(L), tests/0D/alMelting_0D\n"
         "    type            eMelting;\n    Cv              1059;\n    CvLiquid        1177;\n"
         "    Tm              933.47;\n    Lm              3.966e5;\n    dTm             2;\n    Hf              0;")

    # Al burning: Beckstead with ignition temperature
    pp = os.path.join(c, "phaseProperties")
    with open(pp) as f:
        s = f.read()
    old = s[s.index("        reactionRate    pressureBased;"):s.index("        YoxidiserRef    0.24;\n") + len("        YoxidiserRef    0.24;\n")]
    tign = f"{a.Tign:g}" if a.variant == "full" else "1e5"
    s = s.replace(old, f"""        reactionRate    Beckstead;
        // Beckstead (2005) d^n law; ignition at Tign (1e5: never ignites)
        Tign            {tign};
        // Burns only with gas-phase oxidiser (X_eff of O2, H2O, CO2), scaled
        // by the gas fraction of the surrounding fluid, to Al2O3 + H2/CO
        oxidiserPhase   gas;
        products        true;
        // Ignition tracked per particle cloud (transported ignited fraction)
        ignitedFraction true;
        // Burn down to dilute clouds (the default 1e-4 is reached at 0.5 m)
        alphaMin        1e-6;
""")
    # drag and heat transfer of Al in gas
    s = re.sub(r"(    \(Al in gas\)\n    \{\n        type            )GidaspowErgunWenYu;\n        residualRe      1e-3;\n        alphaCutoff     1e-2;",
               r"\1Loth;\n        alphaCutoff     1e-6;", s, count=1)
    s = s.replace("    solveODE        true;\n    alphaCutoff     1e-2;", "    solveODE        true;\n    alphaCutoff     1e-6;")
    s = re.sub(r"(    \(Al in gas\)\n    \{\n        type            RanzMarshall;\n        NuModel         )Gunn;",
               r"\1Kavanau;", s, count=1)
    s = s.replace("""        type            reactingParticle;
        productSpecies
        (
            (Al 1.0)
        );""", """        type            reactingParticle;
        // Products and consumed oxidisers from the Beckstead rate
        productsFromReactionRate true;""")
    with open(pp, "w") as f:
        f.write(s)
    for key in ("Loth;", "Kavanau;", "Beckstead;", "alphaCutoff     1e-6;", "productsFromReactionRate true;"):
        if key not in s:
            sys.exit(f"phaseProperties: {key} not set")

    # Gas species: add the Al sub-oxides and liquid Al2O3
    pg = os.path.join(c, "physicalProperties.gas")
    with open(pg) as f:
        s = f.read()
    names = re.search(r"\nspecies\n\((.*?)\);", s, re.S).group(1).split()
    names = [n.strip('"') for n in names]
    new = names + [n for n in EXTRA_SPECIES if n not in names]
    s = s.replace(re.search(r"\nspecies\n\((.*?)\);", s, re.S).group(0),
                  "\nspecies\n(\n" + "".join(f"    {n if '(' not in n else chr(34) + n + chr(34)}\n" for n in new) + ");")
    with open(pg, "w") as f:
        f.write(s)
    base = os.path.join(BLAST, "tests/1D_spherical/ChargeCompBAl_1D_spherical_reacting/constant/thermo.compressibleGas")
    subprocess.run([sys.executable, os.path.join(BLAST, "utilities/thermo_fitter/nasa9_species.py"),
                    "--base", base, "--species", *new, "-o", os.path.join(c, "thermo.nasa9Gas")],
                   check=True, capture_output=True)
    subprocess.run([sys.executable, os.path.join(BLAST, "utilities/mechanism/make_gas_mechanism.py"),
                    "--species", *new, "-o", os.path.join(c, "reactions.gas")],
                   check=True, capture_output=True)

    # Products -> gas transfer
    fv = os.path.join(c, "fvModels")
    if a.variant == "inert":
        shutil.move(fv, fv + ".disabled")
    else:
        if a.pThreshold:
            edit(fv, "    pThreshold      1e7;", f"    pThreshold      {a.pThreshold:g};")
            edit(fv, "    pRange          1e6;", f"    pRange          {0.1*a.pThreshold:g};")
        if a.tau:
            edit(fv, "    tau             2e-5;", f"    tau             {a.tau:g};")
    print(f"E2 variant {a.variant}: gas species {len(new)}, Tign {tign}"
          + (f", pThreshold {a.pThreshold:g}" if a.pThreshold else "") + (f", tau {a.tau:g}" if a.tau else ""))


if __name__ == "__main__":
    main()
