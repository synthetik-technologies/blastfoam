# Comparison with CONWEP (Kingery-Bulmash) Results

## Notes

This basic case is designed to provide comparison results with the CONWEP airblast load calculator, which leverages the Simplified Kingery-Bulmash Airblast Equations (https://apps.dtic.mil/dtic/tr/fulltext/u2/a526744.pdf)

A hemispherical surface burst detonation event is considered with a 10kg hemi-spherical TNT charge.  The charge is center detonated, using a linear activation model which based on a predetermined detonation velocity.

The case is modeled using an axisymmetic wedge with a computational domain size of 16x16m.

The time of arrival, maximum overpressure, and maximum impulse are compared between 1 m and 10 m. The time of arrival is determined by the time of the peak pressure.

During/after running the simulation, probe outputs are available in the postProcessing directory.
