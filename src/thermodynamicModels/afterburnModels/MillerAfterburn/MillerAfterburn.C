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
        0.0
    ),
    pName_(dict_.lookupOrDefault("pName", word("p"))),
    p_(mesh_.lookupObject<volScalarField>(pName_)),
    alphaRhoName_
    (
        dict.lookupOrDefault
        (
            "rhoName",
            phaseName == word::null
          ? "rho"
          : IOobject::groupName("alphaRho", phaseName)
        )
    ),
    alphaRhoPhiName_
    (
        dict.lookupOrDefault
        (
            "rhoPhiName",
            phaseName == word::null
          ? "rhoPhi"
          : IOobject::groupName("alphaRhoPhi", phaseName)
        )
    ),

    Q0_("Q0", sqr(dimVelocity), dict_),
    m_(readScalar(dict.lookup("m"))),
    n_(readScalar(dict.lookup("n"))),
    a_("a", pow(dimPressure, -n_)/dimTime, dict_),
    pMin_("pMin", dimPressure, dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::afterburnModels::MillerAfterburn::~MillerAfterburn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::afterburnModels::MillerAfterburn::setODEFields
(
    const label nSteps,
    const labelList& oldIs,
    const label& nOld,
    const labelList& deltaIs,
    const label nDelta
)
{
    afterburnModel::setODEFields(nSteps, oldIs, nOld, deltaIs, nDelta);
    cOld_.resize(nOld_);
    deltaC_.resize(nDelta_);
    deltaAlphaRhoC_.resize(nDelta_);
}


void Foam::afterburnModels::MillerAfterburn::clearODEFields()
{
    cOld_.clear();
    cOld_.resize(nOld_);

    ddtC_.clear();

    deltaC_.clear();
    deltaC_.resize(nDelta_);

    deltaAlphaRhoC_.clear();
    deltaAlphaRhoC_.resize(nDelta_);
}


void Foam::afterburnModels::MillerAfterburn::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    const volScalarField& alphaRho
    (
        c_.mesh().lookupObject<volScalarField>(alphaRhoName_)
    );
    const surfaceScalarField& alphaRhoPhi
    (
        c_.mesh().lookupObject<surfaceScalarField>(alphaRhoPhiName_)
    );

    if (oldIs_[stepi - 1] != -1)
    {
        cOld_.set
        (
            oldIs_[stepi - 1],
            new volScalarField(c_)
        );
    }
    volScalarField cOld(ai[stepi - 1]*c_);
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
           cOld += ai[fi]*cOld_[fi];
        }
    }

    tmp<volScalarField> p(p_*pos(p_ - pMin_));
    p.ref().max(small);
    volScalarField deltaC
    (
        a_*pow(max(1.0 - c_, 0.0), m_)*pow(p, n_)
    );
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaC_.set
        (
            deltaIs_[stepi - 1], new volScalarField(deltaC)
        );
    }

    deltaC *= bi[stepi - 1];
    scalar f = bi[stepi - 1];
    for (label i = 0; i < stepi - 1; i++)
    {
        f += bi[i];
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaC += bi[fi]*deltaC_[fi];
        }
    }
    dimensionedScalar dT = alphaRho.time().deltaT();
    c_ = cOld + dT*deltaC;
    c_.min(1);
    c_.max(0);
    c_.correctBoundaryConditions();

    if (!ddtC_.valid())
    {
        ddtC_ = tmp<volScalarField>
        (
            new volScalarField(Foam::max(c_ - cOld, 0.0)/(dT*f))
        );
    }
    else
    {
        ddtC_.ref() = Foam::max(c_ - cOld, 0.0)/(dT*f);
    }

    volScalarField deltaAlphaRhoC(fvc::div(alphaRhoPhi, c_));
    if (deltaIs_[stepi - 1] != -1)
    {
        deltaAlphaRhoC_.set
        (
            deltaIs_[stepi - 1],
            new volScalarField(deltaAlphaRhoC)
        );
    }

    deltaAlphaRhoC *= bi[stepi - 1];
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            deltaAlphaRhoC.ref() +=
                bi[fi]*deltaAlphaRhoC_[fi];
        }
    }
    c_ =
        (
            cOld*alphaRho.oldTime()
            + dT*(ddtC_()*f*alphaRho - deltaAlphaRhoC)
        )/max(alphaRho, dimensionedScalar(dimDensity, 1e-10));
    c_.min(1);
    c_.max(0);
    c_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField>
Foam::afterburnModels::MillerAfterburn::ESource() const
{
    return ddtC_()*Q0_;
}

// ************************************************************************* //
