/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "fluxSchemeBase.H"
#include "MUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxSchemeBase, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemeBase::fluxSchemeBase(const fvMesh& mesh, const word& phaseName)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName("fluxScheme", phaseName),
            mesh.time().timeName(),
            mesh
        )
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemeBase::~fluxSchemeBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::fluxSchemeBase::interpolate
(
    const volScalarField& f,
    const word& fName
) const
{
    autoPtr<MUSCLReconstructionScheme<scalar>> fLimiter
    (
        MUSCLReconstructionScheme<scalar>::New(f, fName)
    );

    tmp<surfaceScalarField> fOwnTmp(fLimiter->interpolateOwn());
    tmp<surfaceScalarField> fNeiTmp(fLimiter->interpolateNei());

    const surfaceScalarField& fOwn = fOwnTmp();
    const surfaceScalarField& fNei = fNeiTmp();

    tmp<surfaceScalarField> tmpff
    (
        surfaceScalarField::New
        (
            fName + "f",
            mesh_,
            dimensioned<scalar>("0", f.dimensions(), Zero)
        )
    );
    surfaceScalarField& ff = tmpff.ref();
    const bool isDensity = (f.dimensions() == dimDensity);

    forAll(fOwn, facei)
    {
        ff[facei] =
            interpolate(fOwn[facei], fNei[facei], isDensity, facei);
    }

    forAll(ff.boundaryField(), patchi)
    {
        scalarField& pff = ff.boundaryFieldRef()[patchi];
        const scalarField& pfOwn = fOwn.boundaryField()[patchi];
        const scalarField& pfNei = fNei.boundaryField()[patchi];
        forAll(pff, facei)
        {
            pff[facei] =
                interpolate
                (
                    pfOwn[facei],
                    pfNei[facei],
                    isDensity,
                    facei, patchi
                );
        }
    }
    return tmpff;
}


Foam::tmp<Foam::surfaceVectorField> Foam::fluxSchemeBase::interpolate
(
    const volVectorField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


Foam::tmp<Foam::surfaceSymmTensorField> Foam::fluxSchemeBase::interpolate
(
    const volSymmTensorField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


Foam::tmp<Foam::surfaceSphericalTensorField>
Foam::fluxSchemeBase::interpolate
(
    const volSphericalTensorField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


Foam::tmp<Foam::surfaceTensorField> Foam::fluxSchemeBase::interpolate
(
    const volTensorField& f,
    const word& fName
) const
{
    return interpolateField(f, fName);
}


bool Foam::fluxSchemeBase::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //
