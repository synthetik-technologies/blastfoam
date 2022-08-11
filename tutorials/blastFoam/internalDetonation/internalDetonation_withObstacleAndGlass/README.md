# Internal detonation with obstacle and glass tutorial

## Notes

This case is used to illustrate how to both setup and use burst patches. The burst patches allow internal faces to act in one way prior to a condition being met, and internal faces after. In this case they are treated as walls prior to the overpressure or impulse reaching their burst condition (4e7 Pa and 1e4 Pa.s) respectively. It should be noted that the first two window burst due to the overpressure condition, while the final window bursts do to impulse. This can be seen due to the lag in time between when the primary pressure wave hits the final window and when it actually bursts.

The mesh is created using the following procedure

<ol>
    <li> Base mesh is created with <strong>blockMesh</strong></li>
    <li> Regions where windows will be placed are refined with <strong>setRefinedFields</strong> </li>
    <li> The outlet patch is added using <strong>createPatch</strong></li>
    <li> The windows are created by adding baffles with the <em>burstCyclicAMI</em> patch type with <strong>createBaffles</strong></li>
    <li> Obstructions uses cellSets selected from running <strong>setRefinedFields</strong> and running <strong>subsetMesh internalCells -patch walls</strong></li>
    <li> The charge is the set using <strong>setRefinedFields</strong></li>
    <li> Finally the case can be run with <strong>blastFoam</strong></li>
</ol>


The case takes approximately 5 min to run on a single core desktop