# Reacting Shock tube

## Notes

This is an illustrative case on how to use blastReactingFoam with chemistry. A single methane reaction is used, with two regions: The left containing 77% Nitrogen and 23& Oxygen at STP. The right region is pure methane at 10 atm and 500K. The laminar combustion model is used in combination with and Euler ODE solver for the chemistry reactions.

The calculation took approx. 20 s to run on a single core laptop at using 500 cells.

