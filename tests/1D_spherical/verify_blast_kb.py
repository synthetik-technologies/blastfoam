#!/usr/bin/env python3
"""
verify_blast_kb.py: free-air blast of a 1 kg spherical charge against
Kingery-Bulmash (plan step E1, Gate 3).

    python3 verify_blast_kb.py tnt_blast      # TNT, no equivalence factor
    python3 verify_blast_kb.py compB_blast    # Comp-B, TNT equivalence

The cases (blastFoam, 1D spherical sector, 10 m, centre-initiated
programmed burn, EuDetCode products EOS, ideal-gas air at 101325 Pa and
288.15 K, no afterburning) record the peak overpressure at every radius (temporalMax) and the
pressure at probes, from which the arrival time and the impulse of the first
positive phase are taken. The reference is kingery_bulmash.py (free-air overpressure; impulse
and arrival time from the hemispherical fits at the equivalent free-air
mass). Comp-B uses the TNT equivalence of UFC 3-340-02 (Composition B:
1.11 for pressure and arrival time, 0.98 for impulse).

Checks (regression bounds for the default 4 mm mesh; the Gate 3 target of
5% in peak overpressure is reported, not enforced):
1. Solver finished; overpressure decays monotonically with radius.
2. Arrival time within 3% of KB for Z >= 3.
3. Peak overpressure within 20% of KB for 1 <= Z <= 8.
4. Positive impulse within 25% of KB for 3 <= Z <= 8.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kingery_bulmash as kb  # noqa: E402

P0 = 101325.0
W = 1.0
EQUIVALENCE = {
    "tnt_blast": (1.0, 1.0),
    "compB_blast": (1.11, 0.98),
}
PROBES = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)


def last_time_dir(path):
    return sorted(os.listdir(path), key=float)[-1]


def load(case):
    pp = os.path.join(case, "postProcessing")
    t = last_time_dir(os.path.join(pp, "radialPeakPressure"))
    pm = np.loadtxt(os.path.join(pp, "radialPeakPressure", t, "p_temporalMax.dat"))
    pr = np.loadtxt(os.path.join(pp, "pressureProbes", "0", "p"))
    return pm, pr


def first_positive_impulse(t, dp):
    """Impulse of the first positive phase (from arrival to the first
    return below ambient). temporalIntegral with positiveOnly also adds the
    slight overpressure (about 100 Pa) that remains near the charge
    afterwards (about 2 Pa s by 25 ms)."""
    a = np.argmax(dp > 0.02*P0)
    neg = np.nonzero(dp[a:] < 0)[0]
    e = a + neg[0] if len(neg) else len(dp)
    return np.trapezoid(dp[a:e], t[a:e])


def at(table, R):
    return table[np.argmin(np.abs(table[:, 1] - R)), 5]


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "tnt_blast"
    eq_p, eq_i = EQUIVALENCE[name]
    case = os.path.join(HERE, name)
    log = os.path.join(case, "log.blastFoam")
    if not os.path.isfile(log) or "\nEnd" not in open(log).read():
        print(f"FAIL: {name}: blastFoam did not finish")
        return 1
    pm, pr = load(case)

    ok = True
    rows = []
    print(f"{name}: 1 kg, TNT equivalence {eq_p} (pressure, arrival) / {eq_i} (impulse)")
    print(f"{'R [m]':>6} {'Z':>5} | {'dp [kPa]':>9} {'KB':>9} {'err':>7} | "
          f"{'I [Pa s]':>9} {'KB':>8} {'err':>7} | {'t_a [ms]':>8} {'KB':>8} {'err':>7}")
    for j, R in enumerate(PROBES):
        Z = R/W**(1.0/3.0)
        p_kb, _, t_kb = kb.free_air(R, eq_p*W)
        _, i_kb, _ = kb.free_air(R, eq_i*W)
        dp = at(pm, R) - P0
        col = pr[:, j + 1] - P0
        arrived = col > 0.02*P0
        t_a = pr[np.argmax(arrived), 0] if arrived.any() else float("nan")
        imp = first_positive_impulse(pr[:, 0], col)
        e_p, e_i, e_t = dp/p_kb - 1, imp/i_kb - 1, t_a/t_kb - 1
        rows.append((R, Z, dp, e_p, e_i, e_t))
        print(f"{R:6.2f} {Z:5.2f} | {dp/1e3:9.1f} {p_kb/1e3:9.1f} {100*e_p:+6.1f}% | "
              f"{imp:9.1f} {i_kb:8.1f} {100*e_i:+6.1f}% | {t_a*1e3:8.3f} {t_kb*1e3:8.3f} {100*e_t:+6.1f}%")

    dps = [r[2] for r in rows]
    if not all(a > b for a, b in zip(dps, dps[1:])):
        print("FAIL: overpressure does not decay monotonically")
        ok = False
    for R, Z, dp, e_p, e_i, e_t in rows:
        if Z >= 3 and abs(e_t) > 0.03:
            print(f"FAIL: arrival time at {R} m differs by {100*e_t:+.1f}% (> 3%)")
            ok = False
        if 1 <= Z <= 8 and abs(e_p) > 0.20:
            print(f"FAIL: peak overpressure at {R} m differs by {100*e_p:+.1f}% (> 20%)")
            ok = False
        if 3 <= Z <= 8 and abs(e_i) > 0.25:
            print(f"FAIL: impulse at {R} m differs by {100*e_i:+.1f}% (> 25%)")
            ok = False

    gate = max(abs(e_p) for R, Z, dp, e_p, e_i, e_t in rows if Z >= 1)
    print(f"Gate 3 (peak overpressure within 5% for Z >= 1): max deviation {100*gate:.1f}% "
          f"({'met' if gate <= 0.05 else 'not met'})")
    print("SUCCESS: blast within the regression bounds of Kingery-Bulmash" if ok else "FAILURE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
