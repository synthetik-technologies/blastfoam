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

#include "activationModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(activationModel, 0);
    defineRunTimeSelectionTable(activationModel, dictionary);
}

// * * * * * * * * * * * * * detonationPoint Functions * * * * * * * * * * * //

Foam::activationModel::detonationPoint::detonationPoint
(
    const volScalarField& alpha,
    const vector& pt,
    const scalar& delay,
    const scalar& radius
)
:
    vector(pt),
    activated_(false),
    delay_(delay),
    radius_(radius)
{
    const fvMesh& mesh = alpha.mesh();

    label nCells = 0;
    scalar sumAlpha = 0.0;
    if (radius > small)
    {
        forAll(mesh.C(), celli)
        {
            if (mag(mesh.C()[celli] - pt) < radius)
            {
                sumAlpha += alpha[celli];
                nCells++;
            }
        }
    }
    else
    {
        label celli = mesh.findCell(pt);
        if (celli >= 0)
        {
            sumAlpha += alpha[celli];
            nCells++;
        }
    }
    if (returnReduce(nCells, sumOp<label>()) == 0)
    {
        WarningInFunction
            << "No cells will be activated using the "
            << "detonation point " << pt << endl;
    }
    if (returnReduce(sumAlpha, sumOp<scalar>()) < small)
    {
        WarningInFunction
            << "No mass was found in the region specified for activation "
            << " for detonation point " << pt << endl;
    }
}

void Foam::activationModel::detonationPoint::setActivated
(
    volScalarField& lambda,
    const bool update
) const
{
    const scalar& t = lambda.time().value();
    if (activated_ || t < delay_)
    {
        return;
    }
    Info<<"activating point " << *this << endl;

    const fvMesh& mesh = lambda.mesh();
    label nCells = 0;
    if (radius_ > small)
    {
        forAll(mesh.C(), celli)
        {
            if (mag(mesh.C()[celli] - *this) < radius_)
            {
                lambda[celli] = 1.0;
                nCells++;
            }
        }
    }
    else
    {
        label celli = mesh.findCell(*this);
        if (celli >= 0)
        {
            lambda[celli] = 1.0;
            nCells++;
        }
    }
    if (returnReduce(nCells, sumOp<label>()) == 0)
    {
        WarningInFunction
            << "No cells were activated using the "
            << "detonation point " << *this << endl;
    }
    if (update)
    {
        activated_ = true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModel::activationModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    integrationSystem
    (
        IOobject::groupName("activationModel", phaseName),
        mesh
    ),
    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0),
        wordList(mesh.boundaryMesh().size(), "zeroGradient")
    ),
    detonationPoints_(0),
    e0_
    (
        dict.found("E0")
      ? dimensionedScalar("E0", dimPressure, dict)
       /dimensionedScalar
        (
            "rho0",
            dimDensity,
            dict.parent().subDict("products").subDict("equationOfState")
        )
      : dimensionedScalar("e0", dimEnergy/dimMass, dict)
    ),
    lambdaExp_(dict.lookupOrDefault("lambdaExp", 1.0)),
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
    alphaRho_(mesh.lookupObject<volScalarField>(alphaRhoName_)),
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
    alphaRhoPhi_(mesh_.lookupObject<surfaceScalarField>(alphaRhoPhiName_)),
    maxDLambda_(dict.lookupOrDefault("maxDLambda", 1.0))
{
    const volScalarField& alpha
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phaseName)
        )
    );

    forAll(detonationPoints_, pti)
    {
        label celli = mesh.findCell(detonationPoints_[pti]);
        if (returnReduce(celli, maxOp<label>()) < 0)
        {
            WarningInFunction
                << "Detonation point at " << detonationPoints_[pti]
                << " is was not found in the mesh. "
                << endl;
        }
        else if (celli >= 0)
        {
            if (alpha[celli] < small)
            {
                WarningInFunction
                    << "There is no mass for phase " << phaseName
                    << " at " << detonationPoints_[pti] << endl;
            }
        }
    }

    Switch useCOM(dict.lookupOrDefault("useCOM", false));
    List<vector> points
    (
        useCOM
      ? List<vector>(1, this->centerOfMass(mesh, alpha))
      : (
            this->needDetonationPoints()
          ? dict.lookup<List<vector>>("points")
          : dict.lookupOrDefault("points", List<vector>(0))
        )
    );

    scalarList delays
    (
        dict.lookupOrDefault
        (
            "delays",
            scalarList(points.size(), 0.0)
        )
    );
    List<scalar> radii
    (
        dict.found("radii")
      ? dict.lookup<List<scalar>>("radii")
      : (
            dict.found("radius")
          ? List<scalar>(points.size(), dict.lookup<scalar>("radius"))
          : List<scalar>(points.size(), 0.0)
        )
    );

    if (points.size())
    {
        if (points.size() != radii.size() || points.size() != delays.size())
        {
            FatalErrorInFunction
                << "User provided detonation points, delays, and radii "
                << "do not match in size." << nl
                << "    # of points: " << points.size() << nl
                << "    # of radii: " << radii.size() << nl
                << "    # of delays: " << delays.size() << nl
                << abort(FatalError);
        }

        Info<< "Initiation Points: " << nl
            << "    " << points << nl
            << "Delays: " << nl
            << "    " << delays << nl
            << "Radii: " << nl
            << "    " << radii << endl;
    }

    detonationPoints_.resize(points.size());
    forAll(detonationPoints_, i)
    {
        detonationPoints_.set
        (
            i,
            new detonationPoint
            (
                alpha,
                points[i],
                delays[i],
                radii[i]
            )
        );
        if (mesh.time().value() > delays[i])
        {
            detonationPoints_[i].activated() = true;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModel::~activationModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::activationModel::centerOfMass
(
    const fvMesh& mesh,
    const volScalarField& alpha
) const
{
    scalarField Vtot(mesh.V()*alpha.primitiveField());
    vectorField m1(Vtot*mesh.C().primitiveField());
    scalar V(gSum(Vtot));
    if (V < small)
    {
        FatalErrorInFunction
            << "No mass was found at the center of mass" << endl
            <<abort(FatalError);
    }

    return gSum(m1)/V;
}


void Foam::activationModel::solve()
{

    // Store old value of lambda, old value of alphaRho is stored in the
    // phaseCompressible system
    volScalarField lambdaOld(lambda_);
    this->storeAndBlendOld(lambdaOld, false);

    dimensionedScalar dT(this->mesh_.time().deltaT());
    dimensionedScalar smallRho("small", dimDensity, 1e-10);

    // Calculate delta due to reaction with no advection
    if (!advection())
    {
        lambda_ = lambdaOld;

        //- Correct the lambda field since zero mass will cause "unactivation"
        //  which is not correct for some models
        forAll(this->detonationPoints_, pointi)
        {
            this->detonationPoints_[pointi].setActivated
            (
                lambda_,
                this->finalStep()
            );
        }
        this->correct();

        // Compute the limited change in lambda
        ddtLambda_ = max(lambda_ - lambdaOld, 0.0)/dT;
        volScalarField& ddtLambda = ddtLambda_.ref();

        //- Compute actual delta for the time step knowing the blended
        ddtLambda = this->calcAndStoreDelta(ddtLambda);

        return;
    }
    // Compute change in lambda with advection included
    else
    {
        volScalarField deltaLambda(this->delta());
        deltaLambda.max(0.0);
        this->storeAndBlendDelta(deltaLambda);

        lambda_ = lambdaOld + deltaLambda*dT;

        // Activate points that are delayed
        forAll(detonationPoints_, pointi)
        {
            detonationPoints_[pointi].setActivated
            (
                lambda_,
                this->finalStep()
            );
        }
        this->correct();
        lambda_.min(1);
        lambda_.max(0);
        lambda_.correctBoundaryConditions();

        // Compute the limited change in lambda
        ddtLambda_ = max(lambda_ - lambdaOld, 0.0)/dT;
        volScalarField& ddtLambda = ddtLambda_.ref();
        volScalarField ddtLambdaLimited(ddtLambda);

        //- Compute actual delta for the time step knowing the blended
        ddtLambda = this->calcAndStoreDelta(ddtLambda);

        volScalarField deltaAlphaRhoLambda(fvc::div(alphaRhoPhi_, lambda_));
        this->storeAndBlendDelta(deltaAlphaRhoLambda);

        //- Solve advection
        lambda_ =
            (
                lambdaOld*alphaRho_.prevIter() - dT*deltaAlphaRhoLambda
            )/max(alphaRho_, smallRho)
          + dT*ddtLambdaLimited;

        //- Correct the lambda field since zero mass will cause "unactivation"
        //  which is not correct for some models
        //  Detonation points are not corrected since they should have mass at
        //  the detonation points
        this->correct();

        lambda_.maxMin(0.0, 1.0);
        lambda_.correctBoundaryConditions();
    }
}


Foam::tmp<Foam::volScalarField> Foam::activationModel::ESource() const
{
    return ddtLambda()*e0_;
}


Foam::tmp<Foam::volScalarField> Foam::activationModel::ddtLambda() const
{
    return ddtLambda_();
}

// ************************************************************************* //
