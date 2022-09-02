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
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        0.0,
        "zeroGradient"
    ),
    pScale_(dict.lookupOrDefault("pScale", 1.0)),
    pName_(dict_.lookupOrDefault("pName", word("p"))),
    p_(mesh_.lookupObject<volScalarField>(pName_)),
    alphaRhoPtr_(nullptr),
    alphaRhoPhiPtr_(nullptr),
    Q0_("Q0", sqr(dimVelocity), dict_),
    m_(readScalar(dict.lookup("m"))),
    n_(readScalar(dict.lookup("n"))),
    a_("a", pow(dimPressure, -n_)/dimTime, dict_),
    pMin_("pMin", dimPressure, dict_)
{
    if (dict.found("tUnits"))
    {
        a_.value() *= pow(10.0, -tUnits[dict.lookup<word>("tUnits")]);
    }
    if (dict.found("pUnits") && !dict.found("pScale"))
    {
        pScale_ = pow(10.0, -pUnits[dict.lookup<word>("pUnits")]);
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
}


void Foam::afterburnModels::MillerAfterburn::solve()
{
    // Do not include volume changes
    volScalarField cOld(c_);
    this->storeAndBlendOld(cOld, false);

    tmp<volScalarField> p(p_*pos(p_ - pMin_));
    if (pScale_ != 1.0)
    {
        p.ref() *= pScale_;
    }
    p.ref().max(small);
    volScalarField deltaC
    (
        a_*pow(max(1.0 - c_, 0.0), m_)*pow(p, n_)
    );
    deltaC.max(0.0);
    this->storeAndBlendDelta(deltaC);

    //- Calculate advection
    const volScalarField& alphaRho = alphaRhoPtr_();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhiPtr_();

    volScalarField deltaAlphaRhoC(fvc::div(alphaRhoPhi, c_));
    this->storeAndBlendDelta(deltaAlphaRhoC);


    dimensionedScalar dT = alphaRho.time().deltaT();
    c_ = cOld + dT*deltaC;
    c_.maxMin(0.0, 1.0);
    c_.correctBoundaryConditions();

    // Compute the limited change in c
    ddtC_ = (Foam::max(c_ - cOld, 0.0)/dT);
    volScalarField& ddtC = ddtC_.ref();

    //- Compute actual delta for the time step knowing the blended value
    //  Not limited to 0 since the delta coefficients can be negative
    //  and store
    ddtC = this->calcAndStoreDelta(ddtC);

    //- Final update of c
    c_ =
        (
            cOld*alphaRho.prevIter() - dT*deltaAlphaRhoC
        )/max(alphaRho, dimensionedScalar(dimDensity, 1e-10))
      + dT*ddtC;
    c_.maxMin(0.0, 1.0);
    c_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField>
Foam::afterburnModels::MillerAfterburn::ESource() const
{
    return ddtC_()*Q0_;
}

// ************************************************************************* //
