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
#include "fluxScheme.H"
#include "phaseFluxScheme.H"
#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxSchemeBase, 0);
    bool fluxSchemeBase::needEnergyFlux = false;
}


const Foam::fluxSchemeBase& Foam::fluxSchemeBase::findFluxScheme
(
    const surfaceScalarField& phi
)
{
    const fvMesh& mesh = phi.mesh();
    const word phaseFluxName
    (
        IOobject::groupName(fluxSchemeBase::typeName, phi.group())
    );
    if (mesh.foundObject<fluxScheme>(fluxSchemeBase::typeName))
    {
        return mesh.lookupObject<fluxScheme>(fluxSchemeBase::typeName);
    }
    else if (mesh.foundObject<phaseFluxScheme>(phaseFluxName))
    {
        return mesh.lookupObject<phaseFluxScheme>(phaseFluxName);
    }
    else
    {
        // If the phase name is not the primary phase
        // check all of the included phases
        HashTable<const fluxSchemeBase*> schemes
        (
            mesh.lookupClass<fluxSchemeBase>()
        );
        forAllConstIter
        (
            HashTable<const fluxSchemeBase*>,
            schemes,
            iter
        )
        {
            if (iter()->phases().found(phi.group()))
            {
                return *iter();
            }
        }
    }

    FatalErrorInFunction
        << "Could not determine a fluxScheme to use for " << phi.name() << endl
        << abort(FatalError);
    return mesh.lookupObject<phaseFluxScheme>(phaseFluxName);
}


bool Foam::fluxSchemeBase::foundFluxScheme
(
    const surfaceScalarField& phi
)
{
    const fvMesh& mesh = phi.mesh();
    const word phaseFluxName
    (
        IOobject::groupName(fluxSchemeBase::typeName, phi.group())
    );
    if (mesh.foundObject<fluxScheme>(fluxSchemeBase::typeName))
    {
        return true;
    }
    else if (mesh.foundObject<phaseFluxScheme>(phaseFluxName))
    {
        return true;
    }
    else
    {
        // If the phase name is not the primary phase
        // check all of the included phases
        HashTable<const fluxSchemeBase*> schemes
        (
            mesh.lookupClass<fluxSchemeBase>()
        );
        forAllConstIter
        (
            HashTable<const fluxSchemeBase*>,
            schemes,
            iter
        )
        {
            if (iter()->phases().found(phi.group()))
            {
                return true;
            }
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemeBase::fluxSchemeBase(const surfaceScalarField& phi)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName("fluxScheme", phi.group()),
            phi.mesh().time().timeName(),
            phi.mesh()
        )
    ),
    mesh_(phi.mesh()),
    phi_(phi)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemeBase::~fluxSchemeBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::fluxSchemeBase::upwindFlux() const
{
    tmp<surfaceScalarField> tuFlux
    (
        surfaceScalarField::New
        (
            IOobject::groupName("upwindFlux", phi_.group()),
            phi_.mesh(),
            dimensionedScalar(phi_.dimensions(), Zero)
        )
    );
    surfaceScalarField& uFlux = tuFlux.ref();
     forAll(uFlux , facei)
    {
        uFlux[facei] = this->interpolate(1.0, -1.0, facei);
    }

    surfaceScalarField::Boundary& buFlux = uFlux.boundaryFieldRef();
    forAll(buFlux, patchi)
    {
        scalarField& puFlux = buFlux[patchi];
        forAll(puFlux, facei)
        {
            puFlux[facei] = this->interpolate(1.0, -1.0, facei, patchi);
        }
    }
    return tuFlux;
}


template<>
Foam::tmp<Foam::surfaceScalarField> Foam::fluxSchemeBase::interpolate
(
    const volScalarField& f,
    const word& fName
) const
{
    autoPtr<ReconstructionScheme<scalar>> fLimiter
    (
        ReconstructionScheme<scalar>::New(f, fName)
    );

    tmp<surfaceScalarField> tfOwn;
    tmp<surfaceScalarField> tfNei;
    fLimiter->interpolateOwnNei(tfOwn, tfNei);

    const surfaceScalarField& fOwn = tfOwn();
    const surfaceScalarField& fNei = tfNei();

    tmp<surfaceScalarField> tff
    (
        surfaceScalarField::New
        (
            fName + "f",
            mesh_,
            dimensioned<scalar>("0", f.dimensions(), Zero)
        )
    );
    surfaceScalarField& ff = tff.ref();

    forAll(fOwn, facei)
    {
        ff[facei] = interpolate(fOwn[facei], fNei[facei], facei);
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
                    pfOwn[facei], pfNei[facei],
                    facei, patchi
                );
        }
    }
    return tff;
}


bool Foam::fluxSchemeBase::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //
