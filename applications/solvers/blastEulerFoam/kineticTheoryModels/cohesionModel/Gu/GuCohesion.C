/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "GuCohesion.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "mathematicalConstants.H"
#include "kineticTheoryModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace cohesionModels
{
    defineTypeNameAndDebug(Gu, 0);

    addToRunTimeSelectionTable
    (
        cohesionModel,
        Gu,
        dictionary
    );

    const Foam::Polynomial<5> Gu::fPhi_
    (
        scalarList
        (
            {
                0.0,
                5e-4,
                0.0048,
                -0.0215,
                0.0249
            }
        )
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::cohesionModels::Gu::Gu
(
    const dictionary& dict,
    const kineticTheoryModel& kt
)
:
    cohesionModel(dict, kt),

    Fmax_("Fmax", dimForce, great),

    alphaMinCohesion_("alphaMinCohesion", dimless, 0),
    alphaMax_("alphaMax", dimless, 0),

    aCoh1_("aCoh1", dimless, 0),
    aCoh2_("aCoh2", dimless, 0),
    aCoh3_("aCoh3", dimless, 0),
    ae_("ae", dimless, 0),
    aw_("aw", dimless, 0),

    tauYield_
    (
        IOobject
        (
            IOobject::groupName("Gu:tauYield", kt.phase().name()),
            kt.phase().mesh().time().timeName(),
            kt.phase().mesh()
        ),
        kt.phase().mesh(),
        dimensionedScalar("0", dimPressure, 0.0),
        extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    W_
    (
        IOobject
        (
            IOobject::groupName("Gu:W", kt.phase().name()),
            kt.phase().mesh().time().timeName(),
            kt.phase().mesh()
        ),
        kt.phase().mesh(),
        dimensionedScalar("0", dimless, 0.0),
        extrapolatedCalculatedFvPatchScalarField::typeName
    )
{
    const dictionary& coeffDict = dict.optionalSubDict(typeName + "Coeffs");
    Fmax_.read(coeffDict);
    alphaMinCohesion_.read(coeffDict);
    alphaMax_.read(coeffDict);
    aCoh1_.read(coeffDict);
    aCoh2_.read(coeffDict);
    aCoh3_.read(coeffDict);
    ae_.read(coeffDict);
    aw_.read(coeffDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::cohesionModels::Gu::~Gu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModels::cohesionModels::Gu::update()
{
    const volScalarField& alpha = kt_.phase();
    const volScalarField& rho = kt_.phase().rho();
    tmp<volScalarField> tTheta = kt_.Theta();
    const volScalarField& Theta = tTheta();
    tmp<volScalarField> td = kt_.phase().d();
    const volScalarField& d = td();
    const scalar e = kt_.es();

    forAll(tauYield_, celli)
    {
        tauYield_[celli] =
            Fmax_.value()*fPhi_.value(max(alpha[celli], small))
           /sqr(d[celli])
           /stabilise(alphaMax_.value() - alpha[celli], 1e-6);
        tauYield_[celli] = max(tauYield_[celli], small);

        const scalar x = aw_.value()*rho[celli]*Theta[celli]/tauYield_[celli];
        const scalar y = sqrt(e*x + 2.0);

        scalar w1 = 2.0*log(1.0 + 0.8842*y);
        W_[celli] =
            (w1 - log(1.0 + 0.9294*log(1.0 + 0.510*y)) - 1.213)
           /(1.0 + w1);
    }

    tauYield_.correctBoundaryConditions();
    W_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::Gu::nu() const
{
    const volScalarField& rho = kt_.phase().rho();
    volSymmTensorField S(dev(symm(fvc::grad(kt_.phase().U()))));
    volScalarField D(sqrt(2.0*(S && S)));
    D.max(small);

    return
        sqr(tauYield_/rho)
       *W_
       /D
       /(aw_*max(kt_.Theta(), dimensionedScalar(sqr(dimVelocity), 1e-6)));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::Gu::Ps() const
{
    const volScalarField& alpha = kt_.phase();
    tmp<volScalarField> td = kt_.phase().d();
    const volScalarField& d = td();

    tmp<volScalarField> tPsCoh
    (
        volScalarField::New
        (
            "Gu:PsCoh",
            alpha.mesh(),
            dimensionedScalar(dimPressure, 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& PsCoh = tPsCoh.ref();

    forAll(PsCoh, celli)
    {
        if (alpha[celli] > alphaMax_.value())
        {
            PsCoh[celli] = 0.0;
        }
        else
        {
            PsCoh[celli] =
              - aCoh1_.value()*Fmax_.value()*sqr(alpha[celli]/d[celli]);
            if (alpha[celli] > alphaMinCohesion_.value())
            {
                PsCoh[celli] +=
                    aCoh2_.value()*Fmax_.value()
                   *sqr((alpha[celli] - alphaMinCohesion_.value())/d[celli])
                   /(alphaMax_.value() - alpha[celli]);
            }
        }
    }
    PsCoh.correctBoundaryConditions();
    return tPsCoh;

}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::Gu::dPsdAlpha() const
{
    const volScalarField& alpha = kt_.phase();
    tmp<volScalarField> td = kt_.phase().d();
    const volScalarField& d = td();

    tmp<volScalarField> tdPsdAlpha
    (
        volScalarField::New
        (
            "Gu:dPsdAlpha",
            alpha.mesh(),
            dimensionedScalar(dimPressure, 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& dPsdAlpha = tdPsdAlpha.ref();

    forAll(dPsdAlpha, celli)
    {
        if (alpha[celli] > alphaMax_.value())
        {
            dPsdAlpha[celli] = 0.0;
        }
        else
        {
            dPsdAlpha[celli] =
              - 2.0*aCoh1_.value()*Fmax_.value()*alpha[celli]/sqr(d[celli]);
            if (alpha[celli] > alphaMinCohesion_.value())
            {
                scalar f =
                    (alpha[celli] - alphaMinCohesion_.value())
                   /(alphaMax_.value() - alpha[celli]);

                dPsdAlpha[celli] +=
                    aCoh2_.value()*Fmax_.value()/sqr(d[celli])
                   *f*(2.0 + f);
            }
        }
    }
    dPsdAlpha.correctBoundaryConditions();
    return tdPsdAlpha;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::Gu::dPsdTheta() const
{
    return
        volScalarField::New
        (
            "Gu:dPsdTheta",
            kt_.phase().mesh(),
            dimensionedScalar(dimPressure/sqr(dimVelocity), 0.0)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::Gu::dissipationSource
(
    const dimensionedScalar& deltaT
) const
{
    const phaseModel& phase = kt_.phase();
    const volScalarField& rho = phase.rho();
    tmp<volScalarField> tTheta = kt_.Theta();
    const volScalarField& Theta = tTheta();
    tmp<volScalarField> td = kt_.phase().d();
    const volScalarField& d = td();

    tmp<volScalarField> gammaGoeff
    (
        volScalarField::New
        (
            "Gu:dissipationSource",
            aCoh3_*sqr(tauYield_)*W_/(aw_*rho*d)
        )
    );
    tmp<volScalarField> ThetaStar
    (
        volScalarField::New
        (
            "ThetaStar",
            pow
            (
                (
                    phase.alphaRho()*pow(Theta, 1.5)
                  + gammaGoeff*deltaT
                )/(phase.alphaRho() + phase.residualAlphaRho()),
                2.0/3.0
            )
        )
    );
    return 1.5*phase.alphaRho()*(ThetaStar - Theta);

}


// ************************************************************************* //
