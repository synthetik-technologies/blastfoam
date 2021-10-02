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
    rho_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("rho", phaseName)
        )
    ),
    rho0_("rho0", dimDensity, 0.0),
    pMin_("pMin", dimPressure, small)
{
    I_.read(dict);
    if (I_.value() > 0)
    {
        a_ = dict.lookup<scalar>("a");
        b_ = dict.lookup<scalar>("b");
        x_ = dict.lookup<scalar>("x");
        maxLambdaI_ = dict.lookup<scalar>("maxLambdaI");
        rho0_.read(dict.parent().subDict("reactants").subDict("equationOfState"));
        needI_ = true;
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

    // Scale the minimum pressure
    pMin_ *= pScale_;
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

    forAll(R, celli)
    {
        // Remove pressures less than minimum pressure
        scalar p = p_[celli]*pScale_;

        const scalar lambdai = max(lambda_[celli], 0.0);
        const scalar oneMLambda = max(1.0 - lambdai, 0.0);
        if (needI_ && lambdai <= maxLambdaI_)
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

// ************************************************************************* //
