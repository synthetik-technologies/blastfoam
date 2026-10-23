#!/usr/bin/env python3
"""Peak overpressure, arrival time (5% overpressure) and positive impulse at
the pressure probes of the 2 m spherical Comp-B/Al case."""

import os
import re

import numpy as np

p0 = 101325.0
case = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(case, "postProcessing", "pressureProbes", "0", "p")
radii = [float(re.search(r"\(([^ ]+)", line).group(1))
         for line in open(path) if line.startswith("# Probe")]
d = np.loadtxt(path)
t = d[:, 0]

print(f"{'r [m]':>6} {'peak [MPa]':>11} {'t_peak [us]':>12} {'arrival [us]':>13} {'I+ [Pa s]':>10}")
for i, r in enumerate(radii):
    dp = d[:, i + 1] - p0
    hit = dp > 0.05*p0
    if not hit.any():
        print(f"{r:6.2f} {'not reached':>11}")
        continue
    k = np.argmax(hit)
    pos = np.clip(dp[k:], 0, None)
    # positive phase: from arrival to the first return below ambient
    ended = (dp[k:] < 0).any()
    neg = np.argmax(dp[k:] < 0) if ended else len(pos)
    impulse = np.trapezoid(pos[:neg], t[k:k + neg])
    print(f"{r:6.2f} {dp.max()/1e6:11.4g} {t[np.argmax(dp)]*1e6:12.1f} {t[k]*1e6:13.1f} {impulse:10.4g}"
          + ("" if ended else "  (positive phase not finished)"))
