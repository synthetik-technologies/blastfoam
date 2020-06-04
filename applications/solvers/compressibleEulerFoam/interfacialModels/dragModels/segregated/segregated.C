/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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

#include "segregated.H"
#include "phasePair.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(segregated, 0);
    addToRunTimeSelectionTable(dragModel, segregated, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::segregated::segregated
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    m_("m", dimless, dict),
    n_("n", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::segregated::~segregated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::segregated::CdRe
(
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the segregated model."
        << exit(FatalError);

    return pair_.phase1();
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::segregated::K
(
    const label nodei,
    const label nodej
) const
{
    const fvMesh& mesh(pair_.phase1().mesh());

    volScalarField alpha1(pair_.phase1().volumeFraction());
    volScalarField alpha2(pair_.phase2().volumeFraction());

    tmp<volScalarField> trho1(pair_.phase1().rho());
    tmp<volScalarField> trho2(pair_.phase2().rho());

    tmp<volScalarField> tnu1(pair_.phase1().nu());
    tmp<volScalarField> tnu2(pair_.phase2().nu());

    const volScalarField& rho1(trho1());
    const volScalarField& rho2(trho2());

    const volScalarField& nu1(tnu1());
    const volScalarField& nu2(tnu2());

    volScalarField L
    (
        IOobject
        (
            "L",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("L", dimLength, 0),
        zeroGradientFvPatchField<scalar>::typeName
    );

    L.primitiveFieldRef() = cbrt(mesh.V());
    L.correctBoundaryConditions();

    volScalarField I
    (
        alpha1
        /max
        (
            alpha1 + alpha2,
            pair_.phase1().residualAlpha() + pair_.phase2().residualAlpha()
        )
    );

    if (pair_.phase2().nNodes() > 1)
    {
        // Scale so not counted nNodes times
        I *=
            alpha2
           /max
            (
                pair_.phase2(),
                pair_.phase2().residualAlpha()
            );
    }

    volScalarField magGradI
    (
        max
        (
            mag(fvc::grad(I)),
            (pair_.phase1().residualAlpha() + pair_.phase2().residualAlpha())/L
        )
    );

    volScalarField muI
    (
        rho1*nu1*rho2*nu2
       /(rho1*nu1 + rho2*nu2)
    );

    volScalarField muAlphaI
    (
        alpha1*rho1*nu1*alpha2*rho2*nu2
       /(
            max(alpha1, pair_.phase1().residualAlpha())*rho1*nu1
          + max(alpha2, pair_.phase2().residualAlpha())*rho2*nu2
        )
    );

    volScalarField ReI
    (
        pair_.rho()
       *pair_.magUr(nodei, nodej)
       /(magGradI*muI)
    );

    volScalarField lambda(m_*ReI + n_*muAlphaI/muI);

    return lambda*sqr(magGradI)*muI;
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::segregated::Kf
(
    const label nodei,
    const label nodej
) const
{
    return fvc::interpolate(K(nodei, nodej));
}


Foam::scalar Foam::dragModels::segregated::CdRe
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the segregated model."
        << exit(FatalError);

    return pair_.phase1()[celli];
}


Foam::scalar Foam::dragModels::segregated::K
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    NotImplemented;
    return pair_.phase1()[celli];
}



// ************************************************************************* //
