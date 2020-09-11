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
    activationModel(mesh, dict, phaseName),

    I_("I", inv(dimTime), 0.0),
    a_("a", dimless, 0.0),
    b_("b", dimless, 0.0),
    x_("x", dimless, 0.0),
    maxLambdaI_("maxLambdaI", dimless, 1.0),

    G1_("G1", dimPressure, 0.0),
    c_("c", dimless, 0.0),
    d_("d", dimless, 0.0),
    y_("y", dimless, 0.0),
    minLambda1_("minLambda1", dimless, 0.0),
    maxLambda1_("maxLambda1", dimless, 1.0),

    G2_("G2", dimless, 0.0),
    e_("e", dimless, 0.0),
    f_("f", dimless, 0.0),
    z_("z", dimless, 0.0),
    minLambda2_("minLambda2", dimless, 0.0),
    maxLambda2_("maxLambda2", dimless, 1.0),

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
    pMin_("pMin", dimPressure, dict)
{
    I_.read(dict);
    if (I_.value() > 0)
    {
        a_.read(dict);
        b_.read(dict);
        x_.read(dict);
        maxLambdaI_.read(dict);
        rho0_.read(dict.parent().subDict("products").subDict("equationOfState"));
    }

    G1_.read(dict);
    if (G1_.value() > 0)
    {
        c_.read(dict);
        d_.read(dict);
        y_.read(dict);
        minLambda1_.read(dict);
        maxLambda2_.readIfPresent(dict);
        G1_.dimensions().reset(pow(dimPressure, -y_)/dimTime);
    }

    G2_.read(dict);
    if (G2_.value() > 0)
    {
        e_.read(dict);
        f_.read(dict);
        z_.read(dict);
        minLambda2_.read(dict);
        maxLambda2_.readIfPresent(dict);
        G2_.dimensions().reset(pow(dimPressure, -z_)/dimTime);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::pressureBasedActivation::~pressureBasedActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::activationModels::pressureBasedActivation::delta() const
{
    // Remove pressures less than minimum pressure
    tmp<volScalarField> p(p_*pos(p_ - pMin_));
    p.ref().max(small);
    tmp<volScalarField> R
    (
        new volScalarField
        (
            IOobject
            (
                "R",
                p_.mesh().time().timeName(),
                p_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            p_.mesh(),
            dimensionedScalar("0", inv(dimTime), 0.0)
        )
    );
    volScalarField oneMLambda(Foam::max(1.0 - lambda_, 0.0));
    if (I_.value() > 0)
    {
        R.ref() +=
            I_
           *pow(Foam::max(rho_/rho0_ - 1.0 - a_, 0.0), x_)
           *pow(oneMLambda, b_)
           *pos0(maxLambdaI_ - lambda_);
    }
    if (G1_.value() > 0)
    {
        R.ref() +=
            G1_
           *pow(oneMLambda, c_)
           *pow(lambda_, d_)
           *pow(p, y_)
           *pos0(lambda_ - minLambda1_)
           *pos(maxLambda1_ - lambda_);
    }
    if (G2_.value() > 0)
    {
        R.ref() +=
            G2_
           *pow(oneMLambda, e_)
           *pow(lambda_, f_)
           *pow(p, z_)
           *pos0(lambda_ - minLambda2_)
           *pos(maxLambda2_ - lambda_);
    }
    return R;
}

// ************************************************************************* //
