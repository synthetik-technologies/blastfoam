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

#include "pressureBasedActivation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activationModels
{
    defineTypeNameAndDebug(pressureBasedActivation, 0);
    addToRunTimeSelectionTable(activationModel, pressureBasedActivation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModels::pressureBasedActivation::pressureBasedActivation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    activationModel(mesh, dict, phaseName, false),

    pScale_(dict.lookupOrDefault("pScale", 1.0)),

    I_("I", inv(dimTime), 0.0),
    a_(0.0),
    b_(0.0),
    x_(0.0),
    maxLambdaI_(1.0),
    needI_(false),

    G1_("G1", dimPressure, 0.0),
    c_(0.0),
    d_(0.0),
    y_(0.0),
    minLambda1_(-1.0),
    maxLambda1_(1.0),
    needG1_(false),

    G2_("G2", dimless, 0.0),
    e_(0.0),
    f_(0.0),
    z_(0.0),
    minLambda2_(-1.0),
    maxLambda2_(1.0),
    needG2_(false),

    pName_(dict.lookupOrDefault("pName", word("p"))),
    p_(mesh.lookupObject<volScalarField>(pName_)),
    alphaPtr_
    (
        phaseName != word::null
      ? &mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phaseName)
        )
      : nullptr
    ),
    rho_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("rho", phaseName)
        )
    ),
    rho0_("rho0", dimDensity, 0.0),
    pMin_("pMin", dimPressure, small),
    residualAlpha_(dict.lookupOrDefault<scalar>("residualAlpha", 1e-10))
{
    I_.readIfPresent(dict);
    if (I_.value() > 0)
    {
        a_ = dict.lookup<scalar>("a");
        b_ = dict.lookup<scalar>("b");
        x_ = dict.lookup<scalar>("x");
        maxLambdaI_ = dict.lookup<scalar>("maxLambdaI");
        const dictionary& rDict(dict.parent().subDict("reactants"));
        const dictionary& pDict(dict.parent().subDict("products"));
        if (rDict.subDict("equationOfState").found(rho0_.name()))
        {
            rho0_.read(rDict.subDict("equationOfState"));
        }
        else if (pDict.subDict("equationOfState").found(rho0_.name()))
        {
            rho0_.read(pDict.subDict("equationOfState"));
        }
        else
        {
            FatalErrorInFunction
                << "'rho0' was not found if products or reactants" << endl
                << abort(FatalError);
        }
        needI_ = true;
    }
    else if (dict.found("Pcj") && dict.found("vDet"))
    {
        scalar vDet = dict.lookup<scalar>("vDet");
        scalar Pcj = dict.lookup<scalar>("Pcj");
        const dictionary& rDict(dict.parent().subDict("reactants"));
        const dictionary& pDict(dict.parent().subDict("products"));
        if (rDict.subDict("equationOfState").found(rho0_.name()))
        {
            rho0_.read(rDict.subDict("equationOfState"));
        }
        else if (pDict.subDict("equationOfState").found(rho0_.name()))
        {
            rho0_.read(pDict.subDict("equationOfState"));
        }
        else
        {
            FatalErrorInFunction
                << "'rho0' was not found if products or reactants" << endl
                << abort(FatalError);
        }
        I_.value() = 1.0/(1.0/rho0_.value() - Pcj/sqr(rho0_.value()*vDet));
    }
    else
    {
        I_.read(dict);
    }

    G1_.read(dict);
    if (G1_.value() > 0)
    {
        c_ = dict.lookup<scalar>("c");
        d_ = dict.lookup<scalar>("d");
        y_ = dict.lookup<scalar>("y");
        minLambda1_ = dict.lookupOrDefault<scalar>("minLambda1", -1.0);
        maxLambda1_ = dict.lookupOrDefault<scalar>("maxLambda1", 1.0);
        G1_.dimensions().reset(pow(dimPressure, -y_)/dimTime);
        needG1_ = true;
    }

    G2_.read(dict);
    if (G2_.value() > 0)
    {
        e_ = dict.lookup<scalar>("e");
        f_ = dict.lookup<scalar>("f");
        z_ = dict.lookup<scalar>("z");
        minLambda2_ = dict.lookupOrDefault<scalar>("minLambda2", -1.0);
        maxLambda2_ = dict.lookupOrDefault<scalar>("maxLambda2", 1.0);
        G2_.dimensions().reset(pow(dimPressure, -z_)/dimTime);
        needG2_ = true;
    }
    pMin_.readIfPresent(dict);

    // Scale the minimum pressure
    pMin_ *= pScale_;

    if (dict.lookupOrDefault("solveODE", false))
    {
        const dictionary& odeDict = dict.subDict("odeCoeffs");
        solver_ = ODESolver::New(*this, odeDict);
        deltaTDet_.set
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::groupName("pressureBased::deltaT", phaseName),
                    mesh.time().name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    dimTime,
                    odeDict.lookupOrDefault<scalar>
                    (
                        "initialDeltaT",
                        mesh.time().deltaTValue()
                    )
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::pressureBasedActivation::~pressureBasedActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::activationModels::pressureBasedActivation::delta() const
{
    tmp<volScalarField> tR
    (
        volScalarField::New
        (
            IOobject::groupName("pressureBased:R", lambda_.group()),
            p_.mesh(),
            dimensionedScalar("0", inv(dimTime), 0.0)
        )
    );
    volScalarField& R = tR.ref();

    if (solver_.valid())
    {
        const scalar deltaT = mesh().time().deltaTValue();
        scalarField l(nEqns());
        forAll(R, celli)
        {
            if (alphaPtr_.valid() && alphaPtr_()[celli] < residualAlpha_)
            {
                continue;
            }

            const scalar lambdaOld = lambda_[celli];
            if (lambdaOld < 1.0)
            {
                l[0] = lambda_[celli];

                scalar timeLeft = deltaT;
                scalar subDeltaT = min(deltaTDet_()[celli], deltaT);
                while (timeLeft > small)
                {
                    solver_->solve(0, timeLeft, l, celli, subDeltaT);
                    timeLeft -= subDeltaT;
                }
                deltaTDet_()[celli] = subDeltaT;
                R[celli] = (max(min(l[0], 1.0), 0.0) - lambdaOld)/deltaT;
            }
        }
        return tR;
    }

    const volScalarField& alphaRho = alphaRhoPtr_();

    forAll(R, celli)
    {
        // Remove pressures less than minimum pressure
        scalar p = p_[celli]*pScale_;

        const scalar lambdai = max(lambda_[celli], 0.0);
        const scalar oneMLambda = max(1.0 - lambdai, 0.0);
        if (needI_ && alphaRho[celli] > small && lambdai <= maxLambdaI_)
        {
            R[celli] =
                I_.value()
               *pow
                (
                    max(rho_[celli]/rho0_.value() - 1.0 - a_, 0.0),
                    x_
                )
               *pow(oneMLambda, b_);
        }
        if (p > pMin_.value())
        {
            if (needG1_ && lambdai >= minLambda1_ && lambdai < maxLambda1_)
            {
                R[celli] +=
                    G1_.value()
                   *pow(oneMLambda, c_)*pow(lambdai, d_)*pow(p, y_);
            }
            if (needG2_ && lambdai >= minLambda2_ && lambdai < maxLambda2_)
            {
                R[celli] +=
                    G2_.value()
                   *pow(oneMLambda, e_)*pow(lambdai, f_)*pow(p, z_);
            }
        }
        R[celli] = max(R[celli], 0.0);
    }
    return tR;
}

void Foam::activationModels::pressureBasedActivation::derivatives
(
    const scalar t,
    const scalarField& y,
    const label celli,
    scalarField& dfdx
) const
{
    dfdx = 0.0;

    scalar& R = dfdx[0];

    // Remove pressures less than minimum pressure
    scalar p = p_[celli]*pScale_;

    const scalar lambdai = min(max(lambda_[celli], 0.0), 1.0);
    const scalar oneMLambda = 1.0 - lambdai;
    if (needI_ && lambdai <= maxLambdaI_)
    {
        R +=
            I_.value()
           *pow
            (
                max(rho_[celli]/rho0_.value() - 1.0 - a_, 0.0),
                x_
            )
           *pow(oneMLambda, b_);
    }
    if (p > pMin_.value())
    {
        if (needG1_ && lambdai >= minLambda1_ && lambdai < maxLambda1_)
        {
            R += G1_.value()*pow(oneMLambda, c_)*pow(lambdai, d_)*pow(p, y_);
        }
        if (needG2_ && lambdai >= minLambda2_ && lambdai < maxLambda2_)
        {
            R += G2_.value()*pow(oneMLambda, e_)*pow(lambdai, f_)*pow(p, z_);
        }
    }
}


void Foam::activationModels::pressureBasedActivation::jacobian
(
    const scalar t,
    const scalarField& y,
    const label celli,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    dfdx = 0.0;
    dfdy = 0.0;

    scalar& R = dfdx[0];
    scalar& dRdLambda = dfdy(0, 0);

    // Remove pressures less than minimum pressure
    scalar p = p_[celli]*pScale_;

    const scalar lambdai = min(max(lambda_[celli], 0.0), 1.0);
    const scalar oneMLambda = 1.0 - lambdai;
    if (needI_ && alphaRhoPtr_()[celli] > small && lambdai <= maxLambdaI_)
    {
        R =
            I_.value()
           *pow
            (
                max(rho_[celli]/rho0_.value() - 1.0 - a_, 0.0),
                x_
            )
           *pow(oneMLambda, b_);
        if (b_ != 1.0 && oneMLambda > small)
        {
            dRdLambda += b_*R/oneMLambda;
        }
    }
    if (p > pMin_.value())
    {
        if (needG1_ && lambdai >= minLambda1_ && lambdai < maxLambda1_)
        {
            R +=
                G1_.value()
               *pow(oneMLambda, c_)*pow(lambdai, d_)*pow(p, y_);

            if (c_ != 1.0 && oneMLambda > small)
            {
                dRdLambda += c_*R/oneMLambda;
            }
            if (d_ != 1.0 && lambdai > small)
            {
                dRdLambda += d_*R/lambdai;
            }
        }
        if (needG2_ && lambdai >= minLambda2_ && lambdai < maxLambda2_)
        {
            R += G2_.value()*pow(oneMLambda, e_)*pow(lambdai, f_)*pow(p, z_);

            if (e_ != 1.0 && oneMLambda > small)
            {
                dRdLambda += e_*R/oneMLambda;
            }
            if (f_ != 1.0 && lambdai > small)
            {
                dRdLambda += f_*R/lambdai;
            }
        }
    }
}

// ************************************************************************* //
