#!/usr/bin/env python3
"""
compare_e2.py: compare runs of this case set up with setup_e2.py (plan step
E2): peak overpressure and first-positive-phase impulse at the pressure
probes, relative to the first run.

    python3 compare_e2.py inert products full [...]
"""

import os
import re
import sys

import numpy as np

P0 = 101325.0


def probes(case):
    path = os.path.join(case, "postProcessing", "pressureProbes", "0", "p")
    radii = [float(re.search(r"\(([^ ]+)", line).group(1))
             for line in open(path) if line.startswith("# Probe")]
    return radii, np.loadtxt(path)


def first_phase(t, dp):
    if not (dp > 0.05*P0).any():
        return float("nan"), float("nan"), float("nan"), False
    a = np.argmax(dp > 0.05*P0)
    neg = np.nonzero(dp[a:] < 0)[0]
    e = a + neg[0] if len(neg) else len(dp)
    return dp.max(), t[a], np.trapezoid(dp[a:e], t[a:e]), bool(len(neg))


def main():
    cases = sys.argv[1:]
    data = {c: probes(c) for c in cases}
    radii = data[cases[0]][0]
    print(f"{'r [m]':>6} " + " ".join(f"{c:>26s}" for c in cases))
    base = {}
    for j, R in enumerate(radii):
        cells = []
        for c in cases:
            _, d = data[c]
            p, ta, imp, done = first_phase(d[:, 0], d[:, j + 1] - P0)
            if c == cases[0]:
                base[R] = (p, imp)
                cells.append(f"{p/1e3:9.1f} kPa {imp:8.1f} Pa s{'' if done else '*'}")
            else:
                cells.append(f"{100*(p/base[R][0]-1):+8.1f}% {100*(imp/base[R][1]-1):+8.1f}%{'' if done else '*'}  ")
        print(f"{R:6.2f} " + " ".join(f"{x:>26s}" for x in cells))
    print("(first column absolute; others relative to it: peak overpressure, impulse; "
          "* positive phase not finished at the end time)")


if __name__ == "__main__":
    main()
