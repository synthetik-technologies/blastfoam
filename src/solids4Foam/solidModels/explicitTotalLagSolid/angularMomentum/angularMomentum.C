/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "angularMomentum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(angularMomentum, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

angularMomentum::angularMomentum
(
    const fvMesh& vm,
    const dictionary& dict
)
:
    mesh_(vm),
    rho_(vm.lookupObject<volScalarField>("rho"))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
angularMomentum::~angularMomentum()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void angularMomentum::AMconservation
(
    GeometricField<vector, fvPatchField, volMesh>& rhsRhoU,
    GeometricField<vector, fvPatchField, volMesh>& rhsRhoU1,
    const GeometricField<vector, fvPatchField, volMesh>& rhsAm,
    const label RKstage
) const
{
    const scalarField& V = mesh_.V();
    const volVectorField& x(mesh_.lookupObject<volVectorField>("x"));
    const volVectorField& rhoU = mesh_.lookupObject<volVectorField>("rhoU");

    const dimensionedScalar deltaT
    (
        "deltaT",
        dimensionSet(0,0,1,0,0,0,0),
        mesh_.time().deltaTValue()
    );

    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf_x
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            "xAM",
            mesh_,
            dimensioned<vector>("x", x.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh>& xAM = tvf_x.ref();

    tmp<GeometricField<vector, fvPatchField, volMesh> > trhoUAM
    (
        GeometricField<vector, fvPatchField, volMesh>::New
        (
            "rhoUAM",
            mesh_,
            dimensioned<vector>("0", rhoU.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh>& rhoUAM = trhoUAM.ref();

    if (RKstage == 0)
    {
        xAM = x.oldTime();
    }

    else if (RKstage == 1)
    {
        xAM = x.oldTime() + (deltaT/2.0)*(rhoU.oldTime()/rho_);
        rhoUAM = rhoU.oldTime() + (deltaT*rhsRhoU1);
        xAM = xAM + ((deltaT*(rhoUAM/rho_))/2.0);
    }

    tensor K_LL = tensor::zero;
    tensor K_LB = tensor::zero;
    scalar K_BB = 0.0;
    vector R_L = vector::zero;

    forAll(mesh_.cells(), cellID)
    {
        K_LL +=
            V[cellID]
           *((xAM[cellID] & xAM[cellID])*tensor::I - (xAM[cellID]*xAM[cellID]));

        K_LB += V[cellID]
           *tensor(0, -xAM[cellID].z(), xAM[cellID].y(), xAM[cellID].z(), 0,
           -xAM[cellID].x(), -xAM[cellID].y(), xAM[cellID].x(), 0);

        K_BB += -V[cellID];

        R_L +=
            (V[cellID]*rhsAm[cellID])
          + ((V[cellID]*rhsRhoU[cellID]) ^ xAM[cellID]);
    }

    if (Pstream::parRun())
    {
        reduce(K_LL, sumOp<tensor>());
        reduce(K_LB, sumOp<tensor>());
        reduce(K_BB, sumOp<scalar>());
        reduce(R_L, sumOp<vector>());
    }

    tensor LHS = K_LL - ((K_LB & K_LB)/K_BB);
    vector RHS = R_L;

    vector lambda = inv(LHS) & RHS;
    vector beta = (-K_LB & lambda)/K_BB;

    forAll(rhsRhoU, celli)
    {
        rhsRhoU[celli] += (lambda ^ xAM[celli]) + beta;
    }

    volVectorField::Boundary& prhsRhoU(rhsRhoU.boundaryFieldRef());
    forAll(prhsRhoU, patchi)
    {
        prhsRhoU[patchi] +=
            (lambda ^ xAM.boundaryField()[patchi]) + beta;
    }
    if (RKstage == 0)
    {
        rhsRhoU1 = rhsRhoU;
    }

    if (Pstream::parRun())
    {
        rhsRhoU.correctBoundaryConditions();
        rhsRhoU1.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void angularMomentum::printGlobalMomentum
(
    const GeometricField<vector, fvPatchField, volMesh>& rhoU,
    const GeometricField<vector, fvPatchField, volMesh>& x
) const
{

    vector rhoUG = vector::zero;
    vector amG = vector::zero;
    scalar vol = gSum(mesh_.V());

    forAll(mesh_.cells(), cellID)
    {
        rhoUG += rhoU[cellID]*mesh_.V()[cellID];
        amG += mesh_.V()[cellID]*(x[cellID] ^ rhoU[cellID]);
    }

    if (Pstream::parRun())
    {
        reduce(rhoUG, sumOp<vector>());
        reduce(amG, sumOp<vector>());
    }

    Info<< "\nPrinting global momentums ..." << nl
        << "Global linear momentum = " << rhoUG/vol << nl
        << "Global angular momentum = " << amG/vol << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
