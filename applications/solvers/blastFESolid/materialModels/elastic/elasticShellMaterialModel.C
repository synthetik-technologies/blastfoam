/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

#include "elasticShellMaterialModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{
    defineTypeNameAndDebug(elasticShell, 0);
    addNamedToRunTimeSelectionTable
    (
        shellMaterialModel,
        elasticShell,
        linear,
        linearElastic
    );
    addNamedToRunTimeSelectionTable
    (
        shellMaterialModel,
        elasticShell,
        nonLinear,
        neoHookeanElastic
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialModels::elasticShell::elasticShell
(
    const dictionary& dict,
    const feMesh1& mesh,
    const Field<vector>& D,
    const Field<vector>& U,
    const Field<vector>& theta,
    const Field<vector>& omega,
    const GeoType geoType
)
:
    shellMaterialModel(dict, mesh, D, U, theta, omega, geoType),
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
            planeStress_
        );
        lambda_.value() = materialModel::lambda_Enu
        (
            E_.value(),
            nu_.value(),
            planeStress_
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
            planeStress_
        );
        lambda_.value() = materialModel::lambda_Kmu
        (
            K_.value(),
            mu_.value(),
            planeStress_
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

Foam::materialModels::elasticShell::~elasticShell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::materialModels::elasticShell::calcSigma
(
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
    const UIndirectList<vector> D0(this->D_, mesh_.elements()[elemi]);
    const UIndirectList<vector> X0(mesh_.nodes(), mesh_.elements()[elemi]);
    const UIndirectList<vector> theta0(theta_, mesh_.elements()[elemi]);

    const element& el = mesh_.elements()[elemi];
    const scalarRectangularMatrix& dshape(mesh_.dshapes()[elemi][rulei]);
    const tensor& invJ = mesh_.invJs()[elemi][rulei];

    const scalar tH = h_*el.ir()[rulei].z()/2.0;
    const vector n(el.calcOrtho(mesh_.Js()[elemi][rulei]));

    List<vector> D(D0);
    List<vector> X(X0);
    forAll(D, i)
    {
        D[i] += tH*(theta0[i] ^ n);
        X[i] += tH*n;
    }

    const vector v(n*el.ir()[rulei].z()/2.0);

    const tensor gradX(refGradient(X, dshape));

    const FixedList<vector, 3> gs(covariantBasis(gradX, v));
    const FixedList<vector, 3> Gs(contravariantBasis(gs));

    const tensor gradD(gradient(D, dshape, invJ));
//     tensor epsilonBar
//     (
//         symmTensor
//         (
//             gradD.x() & gs[0],
//             0.5*((gradD.x() & gs[1]) + (gradD.y() & gs[0])),
//             0.5*((gradD.x() & gs[2]) + (gradD.z() & gs[0])),
//             gradD.y() & gs[1],
//             0.5*((gradD.y() & gs[2]) + (gradD.z() & gs[1])),
//             gradD.z() & gs[2]
//         )
//     );
//
//     const FixedList<vector, 3> es(localBasis(gs));
//
//     symmTensor epsilon(Zero);
//     for (label i = 0; i < 3; i++)
//     {
//         for (label j = 0; j < 3; j++)
//         {
//             const scalar eij = epsilonBar(i, j);
//             epsilon[symmTensor::XX] +=
//                 eij*(Gs[i] & es[vector::X])*(Gs[j] & es[vector::X]);
//
//             epsilon[symmTensor::XY] +=
//                 eij*(Gs[i] & es[vector::X])*(Gs[j] & es[vector::Y]);
//
//             epsilon[symmTensor::XZ] +=
//                 eij*(Gs[i] & es[vector::X])*(Gs[j] & es[vector::Z]);
//
//             epsilon[symmTensor::YY] +=
//                 eij*(Gs[i] & es[vector::Y])*(Gs[j] & es[vector::Y]);
//
//             epsilon[symmTensor::YZ] +=
//                 eij*(Gs[i] & es[vector::Y])*(Gs[j] & es[vector::Z]);
//
//             epsilon[symmTensor::ZZ] +=
//                 eij*(Gs[i] & es[vector::Z])*(Gs[j] & es[vector::Z]);
//         }
//     }
//
// {
//     const UIndirectList<vector> D(this->D_, mesh_.elements()[elemi]);
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
// Info<<epsilon<<" "<<epsilon0<<endl;
// }


    //- Compute stress
    sigma =
        2.0*mu_.value()*epsilon
      + lambda_.value()*tr(epsilon)*symmTensor::I;

    return vp();
}

Foam::scalar Foam::materialModels::elasticShell::calcPiola
(
    tensor& Piola,
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
//     const UIndirectList<vector> D(this->D_, mesh_.elements()[elemi]);
//     const tensor F
//     (
//         calcF
//         (
//             D,
//             theta
//             mesh_.dshapes()[elemi][rulei],
//             mesh_.invJs()[elemi][rulei]
//         )
//     );
//     const scalar J = det(F);
//     const tensor invFT(inv(F, J).T());
//
//     const scalar mu = mu_.value();
//     const scalar K = K_.value();
//     const scalar lambda = lambda_.value();
//
//     sigma =
//         (
//             mu*dev(pow(J, -2.0/3.0)*symm(F & F.T()))
//           + 0.5*K*(sqr(J) - 1.0)*I
//         )/J;
//     Piola = mu*(F - invFT) + lambda*(J - 1.0)*J*invFT;
    return vp();
}


// ************************************************************************* //

