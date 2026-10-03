# blastFoam Automated Verification Report

**Summary**: 14/14 Passed (0 Failed) | **Total Duration**: 122.82 s

| Tier | Test Case | Status | Duration (s) | Diagnostics / Notes |
|---|---|---|---|---|
| Tier 0 | **EuDetCode C++ Unit Suite** | PASSED | 13.31 | 169/169 tests passed |
| Tier 0 | **Thermodynamic Fitter (fit_ejwlc_nobleabel)** | PASSED | 0.58 | Delta p = 0.000% exact continuity verified |
| Tier 1 | **phaseTransfer_0D (0D Isochoric eJWLC->Gas)** | PASSED | 1.41 | Mass & machine-precision energy conservation verified (>99% transfer) |
| Tier 1 | **Aluminium Particle Kinetics & Condensation (WP 1.2, 1.4, 1.5)** | PASSED | 0.50 | Melting plateau, ignition delay, Beckstead kinetics & Al2O3 latent heat verified |
| Tier 1 | **High-T Air Thermodynamics & Dissociation (WP 1.3)** | PASSED | 0.45 | NASA-9 polynomials and high-T O2/N2 dissociation verified up to 6000 K |
| Tier 1 | **pressurePhaseChange (1D Verification)** | PASSED | 2.82 | Completed successfully with clean execution |
| Tier 2 | **ChargeCompBAl_1D_AlInert (Canonical 1D Baseline)** | PASSED | 47.27 | Completed successfully with clean execution |
| Tier 2 | **shock_dissociating_air (1D High-Mach Shock Tube)** | PASSED | 1.93 | 1D blastFoam shock jump state (Ma=5.63) & real-gas dissociation verified (< 1.7% error) |
| Tier 2 | **Multiphase Interface Advection (WP 2.2)** | PASSED | 10.06 | Contact discontinuity advection verified across HLLC, AUSM+, Kurganov & Tadmor |
| Tier 2 | **gunIB1D Particulate Shock Benchmarks (WP 2.3)** | PASSED | 4.23 | Cross-solver agreement with gunIB1D on Houim 1, 2, 3 benchmarks verified (< 2% error) |
| Tier 2 | **Planar CJ/ZND Detonation (WP 2.4)** | PASSED | 0.46 | Planar detonation front speed and CJ pressure verified against EuDetCode (< 0.1% error) |
| Tier 2 | **Planar ZND Detonation with Phase Transfer** | PASSED | 2.31 | Planar ZND detonation coupled with pressurePhaseTransfer verified (< 1% front speed, active phase conversion) |
| Tier 3 | **spherical_blast (WP 3.1 1D Spherical Sector 1kg C4)** | PASSED | 15.75 | blastFoam spherical impulse & peak overpressure verified against Kingery-Bulmash (8.54% err @ Z=2.0m) |
| Tier 3 | **ChargeCompBAl_1D_spherical (1D Spherical Multiphase CJ Detonation)** | PASSED | 21.73 | Origin initiation & spherical CJ detonation speed verified (< 0.7% error vs EuDetCode) |
