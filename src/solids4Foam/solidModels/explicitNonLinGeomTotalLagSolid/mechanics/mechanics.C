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

#include "mechanics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mechanics, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

mechanics::mechanics
(
    const volTensorField& F,
    const operations& ops
)
:
    MeshObject<fvMesh, MoveableMeshObject, mechanics>(F.mesh()),
    mesh_(F.mesh()),

    ops_(ops),

    N_("N", mesh_.Sf()/mesh_.magSf()),

    n_
    (
        IOobject
        (
            "n",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_.Sf()/mesh_.magSf()
    ),

    stabRhoU_
    (
        IOobject
        (
            "stabRhoU",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor("0", dimVelocity, Zero)
    ),

    stabTraction_
    (
        IOobject
        (
            "stabTraction",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor("0", inv(dimVelocity), Zero)
    ),

    stretch_
    (
        IOobject
        (
            "stretch",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimless, 1.0)
    )

{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

mechanics::~mechanics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool mechanics::movePoints()
{
    N_ = mesh_.Sf()/mesh_.magSf();
    return true;
}

void mechanics::correctN(const volTensorField& F)
{
    surfaceTensorField FcInv(inv(fvc::interpolate(F)));
    n_ = (FcInv.T() & N_)/(mag(FcInv.T() & N_));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mechanics::correct
(
    const GeometricField<scalar, fvPatchField, volMesh>& pWaveSpeed,
    const GeometricField<scalar, fvPatchField, volMesh>& sWaveSpeed,
    const GeometricField<tensor, fvPatchField, volMesh>& F
)
{
    // Spatial normals
    surfaceTensorField FcInv(inv(fvc::interpolate(F)));
    n_ = (FcInv.T() & N_)/(mag(FcInv.T() & N_));

    // Stretch
    volTensorField C(F.T() & F);
    forAll(mesh_.cells(), celli)
    {
        ops_.eigenStructure(C[celli]);
        vector eigVal = ops_.eigenValue();
        stretch_[celli] = sqrt(cmptMin(eigVal));
    }

    if (Pstream::parRun())
    {
        stretch_.correctBoundaryConditions();
    }

    // Stabilisation matrices
    stabRhoU_ =
        fvc::interpolate(pWaveSpeed)*n_*n_
      + fvc::interpolate(sWaveSpeed)*(I - (n_*n_));

    stabTraction_ =
        (n_*n_)/fvc::interpolate(pWaveSpeed)
      + (I - (n_*n_))/fvc::interpolate(sWaveSpeed);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mechanics::printCentroid() const
{
    vector sum = vector::zero;
    scalar vol = gSum(mesh_.V());

    forAll(mesh_.cells(), cell)
    {
        sum += mesh_.C()[cell]*mesh_.V()[cell];
    }

    if (Pstream::parRun())
    {
        reduce(sum, sumOp<vector>());
    }

    Info << "\nCentroid of geometry = " << sum/vol << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
