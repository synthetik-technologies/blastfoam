/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "elasticMaterialModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{
    defineTypeNameAndDebug(elastic, 0);
    addNamedToRunTimeSelectionTable
    (
        materialModel,
        elastic,
        linear,
        linearElastic
    );
    addNamedToRunTimeSelectionTable
    (
        materialModel,
        elastic,
        nonLinear,
        neoHookeanElastic
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialModels::elastic::elastic
(
    const dictionary& dict,
    const feMesh1& mesh,
    const pointVectorField& D,
    const pointVectorField& U,
    const bool planeStress,
    const GeoType geoType
)
:
    materialModel(dict, mesh, D, U, planeStress, geoType),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    K_("K", dimPressure, 0.0),
    mu_("mu", dimPressure, 0.0),
    lambda_("lambda", dimPressure, 0.0)
{
    if (dict.found("E") && dict.found("nu"))
    {
        E_.read(dict);
        nu_.read(dict);
        materialModel::Kmu_Enu
        (
            E_.value(),
            nu_.value(),
            K_.value(),
            mu_.value(),
            planeStress
        );
        lambda_.value() = materialModel::lambda_Enu
        (
            E_.value(),
            nu_.value(),
            planeStress
        );
    }
    else if (dict.found("K") && dict.found("mu"))
    {
        K_.read(dict);
        mu_.read(dict);
        materialModel::Enu_Kmu
        (
            K_.value(),
            mu_.value(),
            E_.value(),
            nu_.value(),
            planeStress
        );
        lambda_.value() = materialModel::lambda_Kmu
        (
            K_.value(),
            mu_.value(),
            planeStress
        );
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Either E and nu or K and mu must be provided" << endl
            << abort(FatalIOError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materialModels::elastic::~elastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::materialModels::elastic::calcSigma
(
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
    const UIndirectList<vector> D(this->D_, mesh_.elements()[elemi]);
    const symmTensor epsilon
    (
        symm
        (
            gradient
            (
                D,
                mesh_.dshapes()[elemi][rulei],
                mesh_.invJs()[elemi][rulei]
            )
        )
    );

    //- Compute stress
    sigma =
        2.0*mu_.value()*epsilon
      + lambda_.value()*tr(epsilon)*symmTensor::I;

    return vp();
}


Foam::scalar Foam::materialModels::elastic::calcPiola
(
    tensor& Piola,
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
    const UIndirectList<vector> D(this->D_, mesh_.elements()[elemi]);
    const tensor F
    (
        calcF
        (
            D,
            mesh_.dshapes()[elemi][rulei],
            mesh_.invJs()[elemi][rulei]
        )
    );
    const scalar J = det(F);
    const tensor invFT(inv(F, J).T());

    const scalar mu = mu_.value();
    const scalar K = K_.value();
    const scalar lambda = lambda_.value();

    sigma =
        (
            mu*dev(pow(J, -2.0/3.0)*symm(F & F.T()))
          + 0.5*K*(sqr(J) - 1.0)*I
        )/J;
    Piola = mu*(F - invFT) + lambda*(J - 1.0)*J*invFT;
    return vp();
}

Foam::scalar Foam::materialModels::elastic::vp
(
    const label elemi,
    const label rulei
) const
{
    return sqrt(this->elasticModulus(elemi, rulei)/rho_.value());
}


Foam::scalar Foam::materialModels::elastic::bulkModulus
(
    const label elemi,
    const label rulei
) const
{
    return K_.value();
}


Foam::scalar Foam::materialModels::elastic::elasticModulus
(
    const label elemi,
    const label rulei
) const
{
    return K_.value() + 4.0/3.0*mu_.value();
}


Foam::scalar Foam::materialModels::elastic::shearModulus
(
    const label elemi,
    const label rulei
) const
{
    return mu_.value();
}

// ************************************************************************* //

