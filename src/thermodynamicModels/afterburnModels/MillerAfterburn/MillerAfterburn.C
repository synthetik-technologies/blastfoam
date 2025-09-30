/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "MillerAfterburn.H"
#include "fluxSchemeBase.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace afterburnModels
{
    defineTypeNameAndDebug(MillerAfterburn, 0);
    addToRunTimeSelectionTable(afterburnModel, MillerAfterburn, dictionary);


    const HashTable<label> MillerAfterburn::pUnits
    (
        {
            Tuple2<word, label>("Pa", 1),
            Tuple2<word, label>("kPa", 3),
            Tuple2<word, label>("bar", 5),
            Tuple2<word, label>("MPa", 6),
            Tuple2<word, label>("Mbar", 11)
        }
    );
    const HashTable<label> MillerAfterburn::tUnits
    (
        {
            Tuple2<word, label>("us", -6),
            Tuple2<word, label>("ms", -3),
            Tuple2<word, label>("s", 1),
        }
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::afterburnModels::MillerAfterburn::MillerAfterburn
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    afterburnModel(mesh, dict, phaseName),
    c_
    (
        IOobject
        (
            IOobject::groupName("c", phaseName),
            mesh_.time().name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        0.0,
        "zeroGradient"
    ),
    pScale_(dict_.lookupOrDefault("pScale", 1.0)),
    pName_(dict_.lookupOrDefault("pName", word("p"))),
    p_(mesh_.lookupObject<volScalarField>(pName_)),
    alphaRhoPtr_(nullptr),
    alphaRhoPhiPtr_(nullptr),
    Q0_("Q0", sqr(dimVelocity), dict_),
    m_(readScalar(dict_.lookup("m"))),
    n_(readScalar(dict_.lookup("n"))),
    a_("a", pow(dimPressure, -n_)/dimTime, dict_),
    pMin_("pMin", dimPressure, dict_)
{
    if (dict_.found("tUnits"))
    {
        a_.value() *= pow(10.0, -tUnits[dict_.lookup<word>("tUnits")]);
    }
    if (dict_.found("pUnits") && !dict_.found("pScale"))
    {
        pScale_ = pow(10.0, -pUnits[dict_.lookup<word>("pUnits")]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::afterburnModels::MillerAfterburn::~MillerAfterburn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::afterburnModels::MillerAfterburn::initializeModels()
{
    word alphaRhoName =
        phaseName_ == word::null
      ? "rho"
      : IOobject::groupName("alphaRho", phaseName_);
    word alphaRhoPhiName =
        phaseName_ == word::null
      ? "rhoPhi"
      : IOobject::groupName("alphaRhoPhi", phaseName_);

    alphaRhoPtr_.set(&c_.mesh().lookupObject<volScalarField>(alphaRhoName));
    alphaRhoPhiPtr_.set
    (
        &c_.mesh().lookupObject<surfaceScalarField>(alphaRhoPhiName)
    );

    alphaRhoPtr_->mesh().addTemporaryObject
    (
        reconstruction::ownName(alphaRhoPtr_->name())
    );
    alphaRhoPtr_->mesh().addTemporaryObject
    (
        reconstruction::neiName(alphaRhoPtr_->name())
    );
}


void Foam::afterburnModels::MillerAfterburn::update()
{
    const volScalarField& alphaRho = alphaRhoPtr_();
    dimensionedScalar dT(this->mesh().time().deltaT());

    alphaRhoCOld_ = alphaRho*c_;

    volScalarField p(p_*pos(p_ - pMin_));
    if (pScale_ != 1.0)
    {
        p *= pScale_;
    }
    p.max(small);

    ddtC_ = a_*pow(max(1.0 - c_, 0.0), m_)*pow(p, n_);
    ddtC_.ref().max(0.0);

    // Calculate the deltas using the current value
    deltaAlphaRhoC_ =
        fvc::div(alphaRhoPhiPtr_(), c_)
      - ddtC_()*alphaRho;
}


void Foam::afterburnModels::MillerAfterburn::solve()
{
    if (!alphaRhoCOld_.valid())
    {
        return;
    }

    const volScalarField& alphaRho = alphaRhoPtr_();
    dimensionedScalar dT(this->mesh().time().deltaT());
    dimensionedScalar smallRho("small", dimDensity, 1e-6);

    this->storeAndBlendOld(alphaRhoCOld_.ref());
    this->storeAndBlendDelta(deltaAlphaRhoC_.ref());

    // volScalarField deltaLambda(ddtLambda_());
    // this->blendDelta(deltaLambda);

    //- Update lambda to include advection and reaction
    //  d(alpha rho lambda)/dt = alpha rho d(lambda)/dt + lambda d(alpha rho)/dt
    c_ =
        (alphaRhoCOld_ - deltaAlphaRhoC_*dT)
       /max(alphaRho, smallRho);
      // + deltaLambda*dT;
    c_.maxMin(0.0, 1.0);
    c_.correctBoundaryConditions();
}


// void Foam::afterburnModels::MillerAfterburn::solve()
// {
//     const volScalarField& alphaRho = alphaRhoPtr_();
//     dimensionedScalar dT(this->mesh().time().deltaT());
//     dimensionedScalar smallAlphaRho("small", dimDensity, 1e-10);
//
//     // Calculate the deltas using the current value
//     const fluxSchemeBase& flux = fluxSchemeBase::findFluxScheme(alphaRhoPhiPtr_());
//     volScalarField deltaAlphaRhoC
//     (
//         fvc::div(flux.flux(c_, alphaRhoPtr_(), flux.phi(), false))
//     );
//     this->storeAndBlendDelta(deltaAlphaRhoC);
//
//     tmp<volScalarField> p(p_*pos(p_ - pMin_));
//     if (pScale_ != 1.0)
//     {
//         p.ref() *= pScale_;
//     }
//     p.ref().max(small);
//     volScalarField deltaC
//     (
//         a_*pow(max(1.0 - c_, 0.0), m_)*pow(p, n_)
//     );
//     deltaC.max(0.0);
//     this->storeAndBlendDelta(deltaC);
//
//     // Do not include volume changes
//     this->storeAndBlendOld(c_, false);
//     volScalarField cOld(c_);
//
//     c_ += deltaC*dT;
//     c_.maxMin(0.0, 1.0);
//     c_.correctBoundaryConditions();
//
//     // Compute the limited change in c
//     ddtC_ = (Foam::max(c_ - cOld, 0.0)/dT);
//     volScalarField& ddtC = ddtC_.ref();
//
//     //- Compute actual delta for the time step knowing the blended value
//     //  Not limited to 0 since the delta coefficients can be negative
//     //  and store
//     ddtC = this->calcAndStoreDelta(ddtC);
//
//     //- Final update of c
//     c_ =
//         cOld*(2.0 - alphaRho/max(alphaRho.prevIter(), smallAlphaRho))
//       + dT*(deltaC - deltaAlphaRhoC/max(alphaRho.prevIter(), smallAlphaRho));
//     c_.maxMin(0.0, 1.0);
//     c_.correctBoundaryConditions();
// }


Foam::tmp<Foam::volScalarField>
Foam::afterburnModels::MillerAfterburn::ESource() const
{
    return ddtC_()*Q0_;
}

// ************************************************************************* //
