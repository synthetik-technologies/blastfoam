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

#include "neoHookeanElasticPlasticMaterialModel.H"
#include "cubicEqn.H"
#include "syncTools.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{
    defineTypeNameAndDebug(neoHookeanElasticPlastic, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialModels::neoHookeanElasticPlastic::neoHookeanElasticPlastic
(
    const dictionary& dict,
    const feMesh1& mesh,
    const pointVectorField& D,
    const pointVectorField& U,
    const bool planeStress,
    const GeoType geoType
)
:
    elasticPlastic(dict, mesh, D, U, planeStress, geoType),
    bBar_(mesh.nIp(), symmTensor::I),
    DbBar_(mesh.nIp(), symmTensor::zero),
    bBarConsistent_(dict.lookupOrDefault<Switch>("bBarConsistent", true))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materialModels::neoHookeanElasticPlastic::~neoHookeanElasticPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::materialModels::neoHookeanElasticPlastic::calcIBar
(
    const symmTensor& devbBar
) const
{
    const scalar detDevbBar(det(devbBar));
    const scalar dddevbBar(devbBar && devbBar);
    const scalar fac1  = 2.0/3.0*dddevbBar;

    scalar alpha1 = 0.0;
    if (fac1 < small)
    {
        alpha1 = 3.0;
    }
    else
    {
        const scalar fac2 = 4.0*(1.0 - detDevbBar)/pow(fac1, 1.5);
        alpha1 = 3.0*sqrt(fac1);
        if (fac2 >= 1.0)
        {
            alpha1 *= cosh(acosh(fac2)/3.0);
        }
        else if (mag(fac2) < 0.999)
        {
            alpha1 *= cos(acos(fac2)/3.0);
        }
        else if (fac2 < 0)
        {
            alpha1 = -1.0;
        }
    }
    return alpha1/3.0;
}


void Foam::materialModels::neoHookeanElasticPlastic::preUpdate()
{
    DepsilonPEq_ = Zero;
    DepsilonP_ = Zero;
    DsigmaY_ = Zero;
    DbBar_ = Zero;
}


void Foam::materialModels::neoHookeanElasticPlastic::postUpdate
(
    const scalar f
)
{
    if (f != 1.0)
    {
        epsilonPEq_ += f*DepsilonPEq_;
        epsilonP_ += f*DepsilonP_;
        sigmaY_ += f*DsigmaY_;
        bBar_ += f*DbBar_;
    }
    else
    {
        epsilonPEq_ += DepsilonPEq_;
        epsilonP_ += DepsilonP_;
        sigmaY_ += DsigmaY_;
        bBar_ += DbBar_;
    }
}


Foam::scalar Foam::materialModels::neoHookeanElasticPlastic::calcSigma
(
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
    const element& elem = mesh_.elements()[elemi];
    const label ipI = mesh_.ipLabels()[elemi][rulei];

    // Element displacment and velocity
    const UIndirectList<vector> D(this->D_, elem);
    const UIndirectList<vector> U(this->U_, elem);

    // Element data
    const scalarRectangularMatrix& dshape = mesh_.dshapes()[elemi][rulei];
    const tensor& invJ = mesh_.invJs()[elemi][rulei];

    // Deformation gradient tensor
    const tensor F(calcF(D, dshape, invJ));

    // Jacobian of deformation gradient tensor
    const scalar J = det(F);

    // Relative deformation gradient tensor using veclocity
    // \nabla (U) dt = \nabla (\Delta D)
    tensor relFBar
    (
        gradient(U, dshape, invJ)*mesh_.mesh().time().deltaTValue()
      + tensor::I
    );
    relFBar /= cbrt(det(relFBar));

    const scalar mu = mu_.value();
    const scalar K = K_.value();

    const symmTensor& bbar0 = bBar_[ipI];
    const scalar epsPEq0 = epsilonPEq_[ipI];
    // const symmTensor& epsP0 = epsilonP_[ipI];
    const scalar sigY0 = sigmaY_[ipI];
    const symmTensor trialbBar(transform(relFBar, bbar0));

    const scalar magTrialbBar = mag(trialbBar);
    const scalar IBar = tr(trialbBar)/3.0;
    const scalar muBar = mu*IBar;

    const symmTensor sTrial(mu*dev(trialbBar));

    const scalar fTrial = mag(sTrial) - sqrt2By3*J*sigY0;

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
        muBar,
        magTrialbBar,
        det(F)
    );

    const symmTensor DepsP(IBar*Dlam*plasticN);

    // Deviatoric stress
    sigma = sTrial - 2.0*mu*DepsP;

    symmTensor Dbbar(sigma/mu);
    if (bBarConsistent_)
    {
        Dbbar += calcIBar(Dbbar)*symmTensor::I;
    }
    else
    {
        Dbbar += IBar*symmTensor::I;
    }
    Dbbar -= bbar0;

    DbBar_[ipI] = Dbbar;
    DepsilonPEq_[ipI] = sqrt2By3*Dlam;
    DsigmaY_[ipI] = sigY - sigY0;
    DepsilonP_[ipI] = DepsP;

    // Add "pressure"
    sigma += (0.5*K*(sqr(J) - 1.0))*symmTensor::I;
    sigma /= J;


    // Calculate scaling factor
    const scalar scale(1.0);//1.0 - (2.0*muBar*Dlam/max(mag(sTrial), small)));
    return sqrt((scale*(4.0/3.0)*mu + K)/this->rho_.value());
}


Foam::scalar Foam::materialModels::neoHookeanElasticPlastic::calcPiola
(
    tensor& Piola,
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
    const element& elem = mesh_.elements()[elemi];
    const label ipI = mesh_.ipLabels()[elemi][rulei];

    // Element displacment and velocity
    const UIndirectList<vector> D(this->D_, elem);
    const UIndirectList<vector> U(this->U_, elem);

    // Element data
    const scalarRectangularMatrix& dshape = mesh_.dshapes()[elemi][rulei];
    const tensor& invJ = mesh_.invJs()[elemi][rulei];

    // Deformation gradient tensor
    const tensor F(calcF(D, dshape, invJ));

    // Jacobian of deformation gradient tensor
    const scalar J = det(F);

    // Relative deformation gradient tensor using veclocity
    // \nabla (U) dt = \nabla (\Delta D)
    const tensor relF// (tensor::I);
    (
        calcF(U, dshape, invJ)*mesh_.mesh().time().deltaTValue()
      + tensor::I
    );
    const tensor relFBar(relF/cbrt(det(relF)));

    const scalar mu = mu_.value();
    const scalar K = K_.value();

    // const symmTensor bbar0(symm(F & F.T()));
    const symmTensor& bbar0 = bBar_[ipI];
    const scalar epsPEq0 = epsilonPEq_[ipI];
    // const symmTensor& epsP0 = epsilonP_[ipI];
    const scalar sigY0 = sigmaY_[ipI];

    symmTensor trialbBar(transform(relFBar, bbar0));

    const scalar magTrialbBar = mag(trialbBar);
    const scalar IBar = tr(trialbBar)/3.0;
    const scalar muBar = mu*IBar;

    const symmTensor sTrial(mu*dev(trialbBar));

    const scalar fTrial = mag(sTrial) - sqrt2By3*J*sigY0;

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
        muBar,
        magTrialbBar,
        det(F)
    );

    const symmTensor DepsP(IBar*Dlam*plasticN);

    // Deviatoric stress
    sigma = sTrial - 2.0*mu*DepsP;

    symmTensor Dbbar(sigma/mu);
    if (bBarConsistent_)
    {
        Dbbar += calcIBar(Dbbar)*symmTensor::I;
    }
    else
    {
        Dbbar += IBar*symmTensor::I;
    }
    Dbbar -= bbar0;

    DbBar_[ipI] = Dbbar;
    DepsilonPEq_[ipI] = sqrt2By3*Dlam;
    DsigmaY_[ipI] = sigY - sigY0;
    DepsilonP_[ipI] = DepsP;

    // Add "pressure"
    sigma += (0.5*K*(sqr(J) - 1.0))*symmTensor::I;
    sigma /= J;

    Piola = (F.inv() & sigma)*J;

    // Calculate scaling factor
    // const scalar scale(1.0 - (2.0*muBar*Dlam/max(mag(sTrial), small)));
    // return sqrt((scale*(4.0/3.0)*mu + K)/this->rho_.value());

//     scalar scale = 1.0 - 2.0*muBar*Dlam/max(mag(sTrial), small);
    tensor C(F.T() & F);

    scalar minEV = great;

    // Coefficients of the characteristic cubic polynomial (a = 1)
    const scalar b =
      - C.xx() - C.yy() - C.zz();
    const scalar c =
        C.xx()*C.yy() + C.xx()*C.zz() + C.yy()*C.zz()
      - C.xy()*C.yx() - C.yz()*C.zy() - C.zx()*C.xz();
    const scalar d =
      - C.xx()*C.yy()*C.zz()
      - C.xy()*C.yz()*C.zx() - C.xz()*C.zy()*C.yx()
      + C.xx()*C.yz()*C.zy() + C.yy()*C.zx()*C.xz() + C.zz()*C.xy()*C.yx();

    // Solve
    Roots<3> roots = cubicEqn(1, b, c, d).roots();

    // Check the root types
    forAll(roots, i)
    {
        switch (roots.type(i))
        {
            case rootType::real:
                minEV = min(minEV, roots[i]);
                break;
            case rootType::complex:
//                 WarningInFunction
//                     << "Complex eigenvalues detected for tensor: " << C
//                     << endl;
//                 lambda[i] = 0;
                break;
            case rootType::posInf:
//                 lambda[i] = vGreat;
                break;
            case rootType::negInf:
//                 lambda[i] = - vGreat;
                break;
            case rootType::nan:
                FatalErrorInFunction
                    << "Eigenvalue calculation failed for tensor: " << C
                    << exit(FatalError);
        }
    }

    return sqrt((4.0*mu/3.0 + K)/rho_.value())/minEV;
}


void Foam::materialModels::neoHookeanElasticPlastic::write() const
{
    writeIpField("bBar", dimless, bBar_);
    writeIpField("epsilonP", dimless, epsilonP_);
    writeIpField("epsilonPEq", dimless, epsilonPEq_);
    writeIpField("sigmaY", dimPressure, sigmaY_);
}
// ************************************************************************* //

