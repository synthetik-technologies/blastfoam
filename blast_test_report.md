# blastFoam Automated Verification Report

**Summary**: 3/3 Passed (0 Failed) | **Total Duration**: 77.50 s

| Tier | Test Case | Status | Duration (s) | Diagnostics / Notes |
|---|---|---|---|---|
| Tier 3 | **spherical_blast (WP 3.1 1D Spherical Sector 1kg C4)** | PASSED | 15.95 | blastFoam spherical impulse & peak overpressure verified against Kingery-Bulmash (8.54% err @ Z=2.0m) |
| Tier 3 | **ChargeCompBAl_1D_spherical (1D Spherical Multiphase CJ Detonation)** | PASSED | 30.02 | Origin initiation & spherical CJ detonation speed verified (< 0.7% error vs EuDetCode) |
| Tier 3 | **ChargeCompBAl_1D_spherical_reacting (1D Spherical Detonation + Gas Reactions)** | PASSED | 31.53 | Origin initiation, spherical CJ speed & gas-phase reactions verified (< 0.5% error vs EuDetCode) |
