/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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

#include "fluxScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxScheme, 0);
    defineRunTimeSelectionTable(fluxScheme, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxScheme::fluxScheme(const fvMesh& mesh, const word& name)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName("fluxScheme", name),
            mesh.time().timeName(),
            mesh
        )
    ),
    mesh_(mesh),
    dict_(mesh.schemesDict().subDict("fluxSchemes").subDict(name))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxScheme::~fluxScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxScheme::clear()
{
    own_.clear();
    nei_.clear();
    Uf_.clear();
    pf_.clear();
    alphaf_.clear();
}

void Foam::fluxScheme::createSavedFields()
{
    if (own_.valid())
    {
        return;
    }
    own_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("fluxScheme::own", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1.0)
        )
    );
    nei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("fluxScheme::nei", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("-1", dimless, -1.0)
        )
    );
    Uf_ = tmp<surfaceVectorField>
    (
        new surfaceVectorField
        (
            IOobject
            (
                IOobject::groupName("fluxScheme::Uf", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimVelocity, Zero)
        )
    );
    pf_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("fluxScheme::pf", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimPressure, Zero)
        )
    );
    alphaf_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("fluxScheme::alphaf", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimless, Zero)
        )
    );
}

Foam::tmp<Foam::surfaceVectorField> Foam::fluxScheme::Uf() const
{
    if (Uf_.valid())
    {
        return Uf_;
    }
    return tmp<surfaceVectorField>
    (
        new surfaceVectorField
        (
            IOobject
            (
                "fluxScheme::Uf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("0", dimVelocity, Zero)
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::fluxScheme::pf() const
{
    if (pf_.valid())
    {
        return pf_;
    }
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("fluxScheme::pf", this->group()),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", dimPressure, Zero)
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::fluxScheme::alphaf() const
{
    if (alphaf_.valid())
    {
        return alphaf_;
    }
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("fluxScheme::alphaf", this->group()),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", dimless, Zero)
        )
    );
}


void Foam::fluxScheme::update
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaPhi,
    surfaceScalarField& alphaRhoPhi,
    surfaceVectorField& alphaRhoUPhi,
    surfaceScalarField& alphaRhoEPhi
)
{
    createSavedFields();
    const word phaseName(alpha.group());

    surfaceScalarField alphaOwn
    (
        fvc::interpolate(alpha, own_(), scheme("alpha", phaseName))
    );
    surfaceScalarField alphaNei
    (
        fvc::interpolate(alpha, nei_(), scheme("alpha", phaseName))
    );

    surfaceScalarField rhoOwn(fvc::interpolate(rho, own_(), scheme("rho", phaseName)));
    surfaceScalarField rhoNei(fvc::interpolate(rho, nei_(), scheme("rho", phaseName)));

    surfaceVectorField UOwn(fvc::interpolate(U, own_(), scheme("U", phaseName)));
    surfaceVectorField UNei(fvc::interpolate(U, nei_(), scheme("U", phaseName)));

    surfaceScalarField pOwn(fvc::interpolate(p, own_(), scheme("p", phaseName)));
    surfaceScalarField pNei(fvc::interpolate(p, nei_(), scheme("p", phaseName)));

    surfaceScalarField cOwn(fvc::interpolate(c, own_(), scheme("c", phaseName)));
    surfaceScalarField cNei(fvc::interpolate(c, nei_(), scheme("c", phaseName)));

    surfaceScalarField eOwn(fvc::interpolate(e, own_(), scheme("e", phaseName)));
    surfaceScalarField eNei(fvc::interpolate(e, nei_(), scheme("e", phaseName)));

    preUpdate(p);
    forAll(UOwn, facei)
    {

        calculateFluxes
        (
            alphaOwn[facei], alphaNei[facei],
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            alphaPhi[facei],
            alphaRhoPhi[facei],
            alphaRhoUPhi[facei],
            alphaRhoEPhi[facei],
            facei
        );
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {

            calculateFluxes
            (
                alphaOwn.boundaryField()[patchi][facei],
                alphaNei.boundaryField()[patchi][facei],
                rhoOwn.boundaryField()[patchi][facei],
                rhoNei.boundaryField()[patchi][facei],
                UOwn.boundaryField()[patchi][facei],
                UNei.boundaryField()[patchi][facei],
                eOwn.boundaryField()[patchi][facei],
                eNei.boundaryField()[patchi][facei],
                pOwn.boundaryField()[patchi][facei],
                pNei.boundaryField()[patchi][facei],
                cOwn.boundaryField()[patchi][facei],
                cNei.boundaryField()[patchi][facei],
                mesh_.Sf().boundaryField()[patchi][facei],
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoUPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );
        }
    }
    postUpdate();
}

void Foam::fluxScheme::update
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& e,
    const volScalarField& p,
    const volScalarField& c,
    surfaceScalarField& phi,
    surfaceScalarField& alphaRhoPhi,
    surfaceVectorField& alphaRhoUPhi,
    surfaceScalarField& alphaRhoEPhi
)
{
    createSavedFields();

    const word phaseName(alpha.group());

    surfaceScalarField alphaOwn
    (
        fvc::interpolate(alpha, own_(), scheme("alpha", phaseName))
    );
    surfaceScalarField alphaNei
    (
        fvc::interpolate(alpha, nei_(), scheme("alpha", phaseName))
    );

    surfaceScalarField rhoOwn(fvc::interpolate(rho, own_(), scheme("rho", phaseName)));
    surfaceScalarField rhoNei(fvc::interpolate(rho, nei_(), scheme("rho", phaseName)));

    surfaceVectorField UOwn(fvc::interpolate(U, own_(), scheme("U", phaseName)));
    surfaceVectorField UNei(fvc::interpolate(U, nei_(), scheme("U", phaseName)));

    surfaceScalarField pOwn(fvc::interpolate(p, own_(), scheme("p", phaseName)));
    surfaceScalarField pNei(fvc::interpolate(p, nei_(), scheme("p", phaseName)));

    surfaceScalarField cOwn(fvc::interpolate(c, own_(), scheme("c", phaseName)));
    surfaceScalarField cNei(fvc::interpolate(c, nei_(), scheme("c", phaseName)));

    surfaceScalarField eOwn(fvc::interpolate(e, own_(), scheme("e", phaseName)));
    surfaceScalarField eNei(fvc::interpolate(e, nei_(), scheme("e", phaseName)));

    preUpdate(p);
    forAll(UOwn, facei)
    {
        scalar alphaPhi;
        calculateFluxes
        (
            alphaOwn[facei], alphaNei[facei],
            rhoOwn[facei], rhoNei[facei],
            UOwn[facei], UNei[facei],
            eOwn[facei], eNei[facei],
            pOwn[facei], pNei[facei],
            cOwn[facei], cNei[facei],
            mesh_.Sf()[facei],
            phi[facei],
            alphaPhi,
            alphaRhoPhi[facei],
            alphaRhoUPhi[facei],
            alphaRhoEPhi[facei],
            facei
        );
    }

    forAll(U.boundaryField(), patchi)
    {
        forAll(U.boundaryField()[patchi], facei)
        {
            scalar alphaPhi;
            calculateFluxes
            (
                alphaOwn.boundaryField()[patchi][facei],
                alphaNei.boundaryField()[patchi][facei],
                rhoOwn.boundaryField()[patchi][facei],
                rhoNei.boundaryField()[patchi][facei],
                UOwn.boundaryField()[patchi][facei],
                UNei.boundaryField()[patchi][facei],
                eOwn.boundaryField()[patchi][facei],
                eNei.boundaryField()[patchi][facei],
                pOwn.boundaryField()[patchi][facei],
                pNei.boundaryField()[patchi][facei],
                cOwn.boundaryField()[patchi][facei],
                cNei.boundaryField()[patchi][facei],
                mesh_.Sf().boundaryField()[patchi][facei],
                phi.boundaryFieldRef()[patchi][facei],
                alphaPhi,
                alphaRhoPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoUPhi.boundaryFieldRef()[patchi][facei],
                alphaRhoEPhi.boundaryFieldRef()[patchi][facei],
                facei, patchi
            );
        }
    }
    postUpdate();
}


bool Foam::fluxScheme::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //
