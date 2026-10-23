#!/usr/bin/env python3
"""
kingery_bulmash.py: Kingery-Bulmash airblast parameters of TNT.

Simplified Kingery-Bulmash fits of Swisdak (1994) for a hemispherical TNT
surface burst (metric: scaled distance Z = R/W^(1/3) in m/kg^(1/3)):

    value = exp(sum_i c_i (ln Z)^i)

for the incident peak overpressure [kPa], the scaled incident positive
impulse [kPa ms/kg^(1/3)] and the scaled arrival time [ms/kg^(1/3)]. The
coefficients reproduce the upstream validation data
(validation/blastFoam/KingeryBulmash/.../SwisdakData.csv, 10 kg at 1-10 m)
and the examples of IATG 01.80 (2021, Table 5) within 1%.

Free-air (spherical) burst: the incident overpressure uses the Kingery-
Bulmash free-air polynomial (CONWEP form, 0.05 < Z < 40),

    log10(p [kPa]) = sum_{i=0}^{8} A_i U^i,  U = K0 + K1 log10 Z,

with the coefficients as given (in units of 100 kPa, A0 - 2) by
Dos Santos et al. (SIGE 2023, Table I): 935 kPa at Z = 1, 81.6 kPa at Z = 3.
No free-air fits for the impulse and the arrival time were available here;
they are taken from the hemispherical fits for the charge mass W/f(Z) at
which the hemispherical overpressure equals the free-air overpressure at
the same scaled distance (f = 1.32 at Z = 0.3, 1.64-1.67 for Z = 0.75-3,
1.9-2.06 for Z = 8-30; the usual constant reflection factor is 1.8).

References:
    Kingery, C.N., Bulmash, G., 1984. Airblast parameters from TNT spherical
    air burst and hemispherical surface burst. ARBRL-TR-02555.
    Swisdak, M.M., 1994. Simplified Kingery airblast calculations. Minutes
    of the 26th DoD Explosives Safety Seminar.
    Dos Santos et al., 2023. Comparacao entre metodos empiricos para
    obtencao do equivalente TNT. SIGE 2023, ITA (Table I).
"""

import math

import numpy as np

# (Z_max, coefficients) per range; the first range starts at Z_min
_PRESSURE = (0.2, ((2.9, (7.2106, -2.1069, -0.3229, 0.1117, 0.0685)),
                   (23.8, (7.5938, -3.0523, 0.40977, 0.0261, -0.01267)),
                   (198.5, (6.0536, -1.4066))))
_IMPULSE = (0.2, ((0.96, (5.522, 1.117, 0.6, -0.292, -0.087)),
                  (2.38, (5.465, -0.308, -1.464, 1.362, -0.432)),
                  (33.7, (5.2749, -0.4677, -0.2499, 0.0588, -0.00554)),
                  (158.7, (5.9825, -1.062))))
_ARRIVAL = (0.06, ((1.5, (-0.7604, 1.8058, 0.1257, -0.0437, -0.0310, -0.00669)),
                   (40.0, (-0.7137, 1.5732, 0.5561, -0.4213, 0.1054, -0.00929))))


def _fit(Z, table):
    zmin, ranges = table
    if Z < zmin or Z > ranges[-1][0]:
        return float("nan")
    for zmax, c in ranges:
        if Z <= zmax:
            L = math.log(Z)
            return math.exp(sum(ci*L**i for i, ci in enumerate(c)))
    return float("nan")


def hemispherical(R, W):
    """Incident overpressure [Pa], positive impulse [Pa s] and arrival time
    [s] at R [m] from a hemispherical TNT surface burst of W [kg]."""
    w3 = W**(1.0/3.0)
    Z = R/w3
    return (1e3*_fit(Z, _PRESSURE), 1e3*1e-3*_fit(Z, _IMPULSE)*w3, 1e-3*_fit(Z, _ARRIVAL)*w3)


_FREE_K = (-0.214362789151, 1.35034249993)
_FREE_A = (2.611368669, -1.69012801396, 0.00804973591951, 0.336743114941,
           -0.00516226351334, -0.0809228619888, -0.00478507266747,
           0.00793030472242, 0.0007684469735)


def free_air_overpressure(Z):
    """Kingery-Bulmash free-air incident overpressure [Pa] at the scaled
    distance Z [m/kg^(1/3)] (0.05 < Z < 40)."""
    if Z < 0.05 or Z > 40:
        return float("nan")
    U = _FREE_K[0] + _FREE_K[1]*math.log10(Z)
    return 1e3*10**sum(a*U**i for i, a in enumerate(_FREE_A))


def equivalent_factor(Z):
    """Mass ratio f such that the hemispherical overpressure of W/f equals
    the free-air overpressure of W at the same distance (bisection)."""
    target = free_air_overpressure(Z)
    lo, hi = 1.0, 3.0
    for _ in range(80):
        f = 0.5*(lo + hi)
        if hemispherical(Z, 1.0/f)[0] > target:
            lo = f
        else:
            hi = f
    return 0.5*(lo + hi)


def free_air(R, W):
    """Incident overpressure [Pa], positive impulse [Pa s] and arrival time
    [s] at R [m] from a spherical TNT free-air burst of W [kg]."""
    Z = R/W**(1.0/3.0)
    f = equivalent_factor(Z)
    _, i, t = hemispherical(R, W/f)
    return free_air_overpressure(Z), i, t


def free_air_arrays(R, W):
    out = np.array([free_air(r, W) for r in np.atleast_1d(R)])
    return out[:, 0], out[:, 1], out[:, 2]


if __name__ == "__main__":
    for R in (0.5, 1, 2, 4, 8):
        p, i, t = free_air(R, 1.0)
        print(f"1 kg TNT free air, R {R} m: dp {p/1e3:.1f} kPa, i {i:.1f} Pa s, "
              f"t_a {t*1e3:.3f} ms (factor {equivalent_factor(R):.3f})")
