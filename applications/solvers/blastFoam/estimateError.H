// Calculate error across faces based on density gradient
if (!isA<staticFvMesh>(mesh))
{
    vectorField gradRho(fvc::grad(rho)().primitiveField());
    scalarField dL(mesh.V()/fvc::surfaceSum(mesh.magSf()));

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const label nInternalFaces = mesh.nInternalFaces();
    error = 0.0;

    vector solutionD((vector(mesh.solutionD()) + vector::one)/2.0);

    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        label own = owner[facei];
        label nei = neighbour[facei];
        vector dr = mesh.C()[nei] - mesh.C()[own];
        scalar magdr = mag(dr);

        // Ignore error in empty directions
        if (mag(solutionD & (dr/magdr)) > 0.1)
        {
            scalar dRhodr = (rho[nei] - rho[own])/magdr;
            scalar rhoc = (rho[nei] + rho[own])*0.5;
            scalar dl = (dL[own] + dL[nei])*0.5;
            scalar dRhoDotOwn = gradRho[own] & (dr/magdr);
            scalar dRhoDotNei = gradRho[nei] & (-dr/magdr);
            scalar eT =
                max
                (
                    mag(dRhodr - dRhoDotNei)/(0.3*rhoc/dl + mag(dRhoDotNei)),
                    mag(dRhodr - dRhoDotOwn)/(0.3*rhoc/dl + mag(dRhoDotOwn))
                );
            error[own] = max(error[own], eT);
            error[nei] = max(error[nei], eT);
        }
    }
    error.correctBoundaryConditions();
}
