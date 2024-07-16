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
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.pMesh(),
        dimensionedScalar(dimPressure, 0.0)
    ),
    DsigmaY_
    (
        IOobject
        (
            "DsigmaY",
            mesh.time().timeName(),
            mesh.mesh()
        ),
        mesh.pMesh(),
        dimensionedScalar(dimPressure, 0.0)
    ),
    epsilonPEq_
    (
        IOobject
        (
            "epsilonPEq",
            mesh.time().timeName(),
            mesh.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.pMesh(),
        dimensionedScalar(dimless, 0.0)
    ),
    DepsilonPEq_
    (
        IOobject
        (
            "DepsilonPEq",
            mesh.time().timeName(),
            mesh.mesh()
        ),
        mesh.pMesh(),
        dimensionedScalar(dimless, 0.0)
    ),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.pMesh(),
        dimensionedSymmTensor(dimless, Zero)
    ),
    DepsilonP_
    (
        IOobject
        (
            "DepsilonP",
            mesh.time().timeName(),
            mesh.mesh()
        ),
        mesh.pMesh(),
        dimensionedSymmTensor(dimless, Zero)
    )
{
}

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
    const polyMesh& mesh = mesh_.mesh();

    // Sum increments on coupled points
    syncTools::syncPointList
    (
        mesh,
        DepsilonPEq_,
        plusEqOp<scalar>(),
        0.0
    );
    syncTools::syncPointList
    (
        mesh,
        DepsilonP_,
        plusEqOp<symmTensor>(),
        symmTensor::zero
    );
    syncTools::syncPointList
    (
        mesh,
        DsigmaY_,
        plusEqOp<scalar>(),
        0.0
    );

    DepsilonPEq_ /= mesh_.W();
    DepsilonP_ /= mesh_.W();
    DsigmaY_ /= mesh_.W();

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

    //- Make sure coupled points are synced
    mesh_.pushUntransformedData(epsilonPEq_);
    mesh_.pushUntransformedData(epsilonP_);
    mesh_.pushUntransformedData(sigmaY_);

//     epsilonPEq_.correctBoundaryConditions();
//     epsilonP_.correctBoundaryConditions();
//     sigmaY_.correctBoundaryConditions();
}


Foam::scalar Foam::materialModels::elasticPlastic::calcSigma
(
    symmTensor& sigma,
    const label elemi,
    const label rulei
)
{
    const element& elem = mesh_.elements()[elemi];
    const integrationPoint& ip = elem.ir()[rulei];

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

    scalar epsPEq0 = 0.0;
    symmTensor epsP0(Zero);
    scalar sigY0 = 0.0;
    const scalarList& shape = mesh_.shapes()[elemi][rulei];
    {
        const UIndirectList<scalar> epsilonPEq_loc(epsilonPEq_, elem);
        const UIndirectList<symmTensor> epsilonP_loc(epsilonP_, elem);
        const UIndirectList<scalar> sigmaY_loc(sigmaY_, elem);
        forAll(shape, si)
        {
            const scalar s = shape[si];
            epsPEq0 += s*epsilonPEq_loc[si];
            epsP0 += s*epsilonP_loc[si];
            sigY0 += s*sigmaY_loc[si];
        }
    }

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

    const scalar DepsPEq(sqrt2By3*Dlam);
    const scalar DsigY(sigY - sigY0);
    const symmTensor DepsP(Dlam*plasticN);
    if (update_)
    {
        const scalar w = ip.w()*mesh_.Ws()[elemi][rulei];
        UIndirectList<scalar> DepsilonPEq_loc(DepsilonPEq_, elem);
        UIndirectList<symmTensor> DepsilonP_loc(DepsilonP_, elem);
        UIndirectList<scalar> DsigmaY_loc(DsigmaY_, elem);
        forAll(shape, si)
        {
            const scalar sw = shape[si]*w;
            DepsilonPEq_loc[si] += DepsPEq*sw;
            DepsilonP_loc[si] += DepsP*sw;
            DsigmaY_loc[si] += DsigY*sw;
        }
    }

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


// ************************************************************************* //

