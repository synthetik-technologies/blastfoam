/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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
    defineTypeNameAndDebug(activationModel::detonationPoint, 0);
    defineTemplateTypeNameAndDebug
    (
        IOPtrList<activationModel::detonationPoint>,
        0
    );
    defineRunTimeSelectionTable(activationModel, dictionary);
}

// * * * * * * * * * * Static detonationPoint functions  * * * * * * * * * * //

Foam::autoPtr<Foam::activationModel::detonationPoint>
Foam::activationModel::detonationPoint::New(Istream& is)
{
    return autoPtr<detonationPoint>(new detonationPoint(is));
}


// * * * * * * * * * * * * detonationPoint constructor * * * * * * * * * * * //

Foam::activationModel::detonationPoint::detonationPoint
(
    const vector& pt,
    const scalar delay,
    const scalar radius
)
:
    vector(pt),
    activated_(false),
    printedActivated_(false),
    delay_(delay),
    radius_(radius)
{}


Foam::activationModel::detonationPoint::detonationPoint(Istream& is)
:
    vector(is),
    activated_(readBool(is)),
    printedActivated_(false),
    delay_(readScalar(is)),
    radius_(readScalar(is))
{}

// * * * * * * * * * * * * detonationPoint destructor * * * * * * * * * * * //

Foam::activationModel::detonationPoint::~detonationPoint()
{}


// * * * * * * * * * * detonationPoint member functions  * * * * * * * * * * //

Foam::autoPtr<Foam::activationModel::detonationPoint>
Foam::activationModel::detonationPoint::clone() const
{
    return autoPtr<detonationPoint>
    (
        new detonationPoint(*this)
    );
}


bool Foam::activationModel::detonationPoint::check
(
    const volScalarField& alphaRho
)
{
    if (activated_)
    {
        return true;
    }

    const fvMesh& mesh = alphaRho.mesh();
    const vector& pt(*this);

    label nCells = 0;
    scalar sumAlphaRho = 0.0;
    if (radius_ > small)
    {
        forAll(mesh.C(), celli)
        {
            if (mag(mesh.C()[celli] - pt) < radius_)
            {
                sumAlphaRho += alphaRho[celli];
                nCells++;
            }
        }
    }
    else
    {
        label celli = mesh.findCell(pt);
        if (celli >= 0)
        {
            sumAlphaRho += alphaRho[celli];
            nCells++;
        }
    }
    if (returnReduce(nCells, sumOp<label>()) == 0)
    {
        FatalErrorInFunction
            << "No cells will be activated using the "
            << "detonation point " << pt
            << abort(FatalError);
    }
    if (returnReduce(sumAlphaRho, sumOp<scalar>()) < small)
    {
        FatalErrorInFunction
            << "No mass was found in the region specified for "
            << "activation of detonation point " << pt
            << abort(FatalError);
    }
    activated_ = activated_ || mesh.time().value() > delay_;

    return true;
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

    if (!printedActivated_)
    {
        Info<<"activating point " << vector(*this) << endl;
        printedActivated_ = true;
    }
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

    // This is checked again because cells may have been refined or
    // unrefined
    if (returnReduce(nCells, sumOp<label>()) == 0)
    {
        FatalErrorInFunction
            << "No cells were activated using the "
            << "detonation point " << *this
            << abort(FatalError);
    }
    if (update)
    {
        activated_ = true;
    }
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const activationModel::detonationPoint& dp
)
{
    os  << vector(dp) << token::SPACE
        << dp.activated_ << token::SPACE
        << dp.delay_ << token::SPACE
        << dp.radius_ << token::SPACE;
    return os;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModel::activationModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const bool needDetonationPoints
)
:
    timeIntegrationSystem
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
        "zeroGradient"
    ),
    detonationPoints_
    (
        IOobject
        (
            IOobject::groupName("detonationPoints", phaseName),
            mesh.time().timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        readDetonationPoints
        (
            dict,
            mesh.lookupObject<volScalarField>
            (
                IOobject::groupName("alpha", phaseName)
            ),
            needDetonationPoints
        )
    ),
    e0_
    (
        dict.found("E0")
      ? dimensionedScalar("E0", dimPressure, dict)
       /dimensionedScalar
        (
            "rho0",
            dimDensity,
            dict.parent().subDict("products").subDict
            (
                "equationOfState"
            ).found("rho0")
          ? dict.parent().subDict("products").subDict("equationOfState")
          : dict.parent().subDict("reactants").subDict("equationOfState")
        )
      : dimensionedScalar("e0", dimEnergy/dimMass, dict)
    ),
    lambdaExp_(dict.lookupOrDefault("lambdaExp", 1.0)),
    alphaRhoPtr_(nullptr),
    alphaRhoPhiPtr_(nullptr),
    maxDLambda_(dict.lookupOrDefault("maxDLambda", 1.0))
{
    lambda_.storeOldTime();
    if (detonationPoints_.size())
    {
        vectorField points(detonationPoints_.size());
        scalarField delays(detonationPoints_.size());
        scalarField radii(detonationPoints_.size());
        forAll(detonationPoints_, pti)
        {
            points[pti] = detonationPoints_[pti];
            delays[pti] = detonationPoints_[pti].delay();
            radii[pti] = detonationPoints_[pti].radius();
        }
        Info<< "Initiation Points: " << nl
            << "    " << points << nl
            << "Delays: " << nl
            << "    " << delays << nl
            << "Radii: " << nl
            << "    " << radii << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModel::~activationModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PtrList<Foam::activationModel::detonationPoint>
Foam::activationModel::readDetonationPoints
(
    const dictionary& dict,
    const volScalarField& alpha,
    const bool needDetonationPoints
) const
{
    {
        IOobject detPointsHeader
        (
            IOobject::groupName("detonationPoints", alpha.group()),
            alpha.mesh().time().timeName(),
            "uniform",
            alpha.mesh()
        );
        if (detPointsHeader.typeHeaderOk<IOPtrList<detonationPoint>>(true))
        {
            return PtrList<detonationPoint>();
        }
    }

    Switch useCOM(dict.lookupOrDefault("useCOM", false));
    List<vector> points
    (
        (!useCOM || dict.found("points"))
      ? (
            needDetonationPoints
          ? dict.lookup<List<vector>>("points")
          : dict.lookupOrDefault("points", List<vector>(0))
        )
      : List<vector>(1, this->centerOfMass(alpha))
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
    }

    PtrList<detonationPoint> detPoints(points.size());
    forAll(detPoints, i)
    {
        detPoints.set
        (
            i,
            new detonationPoint
            (
                points[i],
                delays[i],
                radii[i]
            )
        );
    }
    return detPoints;
}


void Foam::activationModel::initializeModels()
{
    word phaseName = lambda_.group();
    word alphaRhoName =
        phaseName == word::null
      ? "rho"
      : IOobject::groupName("alphaRho", phaseName);
    word alphaRhoPhiName =
        phaseName == word::null
      ? "rhoPhi"
      : IOobject::groupName("alphaRhoPhi", phaseName);

    if (lambda_.mesh().foundObject<volScalarField>(alphaRhoName))
    {
        alphaRhoPtr_.set
        (
            &lambda_.mesh().lookupObject<volScalarField>
            (
                alphaRhoName
            )
        );

        forAll(detonationPoints_, i)
        {
            detonationPoints_[i].check(alphaRhoPtr_());
        }
    }
    if (lambda_.mesh().foundObject<surfaceScalarField>(alphaRhoPhiName))
    {
        alphaRhoPhiPtr_.set
        (
            &lambda_.mesh().lookupObject<surfaceScalarField>
            (
                alphaRhoPhiName
            )
        );
    }
}

Foam::vector Foam::activationModel::centerOfMass
(
    const volScalarField& alpha
) const
{
    const fvMesh& mesh = alpha.mesh();
    scalarField Vtot(mesh.V()*alpha.primitiveField());
    vectorField m1(Vtot*mesh.C().primitiveField());

    scalar V(gSum(Vtot));
    if (V < small)
    {
        FatalErrorInFunction
            << "No mass was found in the domain"
            << abort(FatalError);
    }

    return gSum(m1)/V;
}


void Foam::activationModel::solve()
{
    dimensionedScalar dT(this->mesh().time().deltaT());
    dimensionedScalar smallRho("small", dimDensity, 1e-10);

    // Calculate the deltas using the current value
    volScalarField deltaAlphaRhoLambda(fvc::div(alphaRhoPhiPtr_(), lambda_));
    this->storeAndBlendDelta(deltaAlphaRhoLambda);

    volScalarField deltaLambda(this->delta());
    deltaLambda.max(0.0);
    this->storeDelta(deltaLambda);

    // Store old value of lambda, old value of alphaRho is stored in the
    // phaseCompressible system
    this->storeAndBlendOld(lambda_, false);
    const volScalarField lambdaOld(lambda_);

    lambda_ += deltaLambda*dT;

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
    lambda_.maxMin(0.0, 1.0);
    lambda_.correctBoundaryConditions();

    // Compute the limited change in lambda
    ddtLambda_ = max(lambda_ - lambdaOld, 0.0)/dT;
    volScalarField& ddtLambda = ddtLambda_.ref();

    //- Solve advection
    lambda_ =
        (
            lambdaOld*alphaRhoPtr_().prevIter()
          - dT*(deltaAlphaRhoLambda - ddtLambda*alphaRhoPtr_())
        )/max(alphaRhoPtr_(), smallRho);

    //- Compute actual delta for the time step knowing the blended
    ddtLambda = this->calcAndStoreDelta(ddtLambda);


    //- Correct the lambda field since zero mass will cause "unactivation"
    //  which is not correct for some models
    //  Detonation points are not corrected since they should have mass at
    //  the detonation points
    this->correct();

    lambda_.maxMin(0.0, 1.0);
    lambda_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField> Foam::activationModel::initESource() const
{
    return volScalarField::New
    (
        "initESource",
        lambda_.mesh(),
        dimensionedScalar("0", e0_.dimensions(), 0.0)
    );
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
