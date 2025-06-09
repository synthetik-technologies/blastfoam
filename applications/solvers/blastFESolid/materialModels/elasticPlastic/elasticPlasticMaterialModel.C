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

#include "elasticPlasticMaterialModel.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{
    defineTypeNameAndDebug(elasticPlastic, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialModels::elasticPlastic::elasticPlastic
(
    const dictionary& dict,
    const feMesh1& mesh,
    const pointVectorField& D,
    const pointVectorField& U,
    const bool planeStress,
    const GeoType geoType
)
:
    elastic(dict, mesh, D, U, planeStress, geoType),
    plasticSolver(dict),
    sigmaY_(mesh_.nIp(), 0.0),
    DsigmaY_(mesh_.nIp(), 0.0),
    epsilonPEq_(mesh_.nIp(), 0.0),
    DepsilonPEq_(mesh_.nIp(), 0.0),
    epsilonP_(mesh_.nIp(), symmTensor::zero),
    DepsilonP_(mesh_.nIp(), symmTensor::zero)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materialModels::elasticPlastic::~elasticPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::materialModels::elasticPlastic::preUpdate()
{
    DepsilonPEq_ = Zero;
    DepsilonP_ = Zero;
    DsigmaY_ = Zero;
}


void Foam::materialModels::elasticPlastic::postUpdate(const scalar f)
{
    if (f != 1.0)
    {
        epsilonPEq_ += f*DepsilonPEq_;
        epsilonP_ += f*DepsilonP_;
        sigmaY_ += f*DsigmaY_;
    }
    else
    {
        epsilonPEq_ += DepsilonPEq_;
        epsilonP_ += DepsilonP_;
        sigmaY_ += DsigmaY_;
    }
}


Foam::scalar Foam::materialModels::elasticPlastic::calcSigma
(
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
    const label ipI = mesh_.ipLabels()[elemi][rulei];

    const UIndirectList<vector> D(this->D_, mesh_.elements()[elemi]);
    const symmTensor eps
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

    const scalar magEps = max(mag(eps), small);

    const scalar mu = mu_.value();
    const scalar K = K_.value();

    const scalar epsPEq0 = epsilonPEq_[ipI];
    const symmTensor& epsP0 = epsilonP_[ipI];
    const scalar sigY0 = sigmaY_[ipI];

    // sTrial = 2*mu(epsilon - dev(epsilonP))
    symmTensor sTrial(2.0*mu*(dev(eps) - dev(epsP0)));

    const scalar fTrial = mag(sTrial) - sqrt2By3*sigY0;

    scalar Dlam = 0.0;
    scalar sigY = sigY0;
    symmTensor plasticN(Zero);
    tmp<yieldStressModel> ysPtr(this->yieldStress(elemi, rulei));
    updatePlasticity
    (
        ysPtr(),
        plasticN,
        Dlam,
        sigY,
        sigY0,
        fTrial,
        sTrial,
        epsPEq0,
        mu,
        magEps,
        1.0
    );

    const symmTensor DepsP(Dlam*plasticN);

    DepsilonPEq_[ipI] = sqrt2By3*Dlam;
    DepsilonP_[ipI] = DepsP;
    DsigmaY_[ipI] = sigY - sigY0;


    sigma = sTrial - 2.0*mu*DepsP + K*tr(eps)*symmTensor::I;

    scalar scale = max(1.0 - 2.0*mu*Dlam/max(mag(sTrial), small), 0.0);
    return sqrt((scale*4.0*mu/3.0 + K)/rho_.value());
}


Foam::scalar Foam::materialModels::elasticPlastic::calcPiola
(
    tensor& Piola,
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
    NotImplemented;
    return 0.0;
}


void Foam::materialModels::elasticPlastic::write() const
{
    writeIpField("epsilonP", dimless, epsilonP_);
    writeIpField("epsilonPEq", dimless, epsilonPEq_);
    writeIpField("sigmaY", dimPressure, sigmaY_);
}
// ************************************************************************* //

