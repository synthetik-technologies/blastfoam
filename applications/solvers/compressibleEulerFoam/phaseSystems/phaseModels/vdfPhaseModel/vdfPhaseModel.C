/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Added population balance model
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "vdfPhaseModel.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vdfPhaseModel::vdfPhaseModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    phaseModel(mesh, dict, phaseName),
    pbeDict_
    (
        IOobject
        (
            "populationBalanceProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    populationBalance_
    (
        populationBalanceModel::New
        (
            name_,
            pbeDict_,
            phiPtr_()
        )
    ),
    quadrature_
    (
        mesh.lookupObjectRef<velocityQuadratureApproximation>
        (
            IOobject::groupName
            (
                "quadratureProperties",
                phaseName
            )
        )
    ),
    computeVariance_(false),
    minD_
    (
        dimensionedScalar::lookupOrDefault
        (
            "minD",
            phaseDict_,
            dimLength,
            1e-6
        )
    ),
    sizeIndex_(quadrature_.nodes()[0].sizeIndex()),
    m0VolumeFraction_(false),
    volumeFractionMoment_(quadrature_.momentOrders()[0].size(), 0),
    sizeMoment_(quadrature_.momentOrders()[0].size(), 0)
{
    const labelList& velocityIndexes =
        quadrature_.nodes()[0].velocityIndexes();
    const labelListList& momentOrders =
        quadrature_.momentOrders();

    forAll(momentOrders, mi)
    {
        forAll(velocityIndexes, cmpt)
        {
            if (momentOrders[mi][velocityIndexes[cmpt]] >= 2)
            {
                computeVariance_ = true;
            }
        }
    }
    if (computeVariance_)
    {
         sigma_ = tmp<volTensorField>
         (
             new volTensorField
            (
                IOobject
                (
                    IOobject::groupName("Sigma", name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedTensor("Sigma", sqr(dimVelocity), Zero)
            )
        );
        Theta_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Theta", name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("Theta", sqr(dimVelocity), 0.0)
            )
        );
    }

    dimensionSet weightDim
    (
        quadrature_.nodes()[0].primaryWeight().dimensions()
    );

    if (weightDim == dimless)
    {
        m0VolumeFraction_ = true;
    }

    if (quadrature_.nodes()[0].sizeIndex() != -1)
    {
        d_.writeOpt() = IOobject::AUTO_WRITE;

        dimensionSet sizeDim
        (
            quadrature_.nodes()[0].primaryAbscissae()[sizeIndex_].dimensions()
        );

        sizeMoment_[sizeIndex_] = 1;

        if (sizeDim == dimLength)
        {
            sizeType_ = length;

            if (weightDim == dimless)
            {
                momentSetType_ = volumeFractionLength;
            }
            else if (weightDim == inv(dimVolume))
            {
                momentSetType_ = numberDensityLength;
                volumeFractionMoment_[sizeIndex_] = 3;
            }
        }
        else if (sizeDim == dimVolume)
        {
            sizeType_ = volume;

            if (weightDim == dimless)
            {
                momentSetType_ = volumeFractionVolume;
            }
            else if (weightDim == inv(dimVolume))
            {
                momentSetType_ = numberDensityVolume;
                volumeFractionMoment_[sizeIndex_] = 1;
            }
        }
        else if (sizeDim == dimMass)
        {
            sizeType_ = mass;

            if (weightDim == dimless)
            {
                momentSetType_ = volumeFractionMass;
            }
            else if (weightDim == inv(dimVolume))
            {
                momentSetType_ = numberDensityMass;
                volumeFractionMoment_[sizeIndex_] = 1;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unknown weight and abscissae dimension combination." << nl
                << "Weight can use dimensions of volume fraction or" << nl
                << "number density and size abscissa can use dimensions" << nl
                << "of mass, length or volume."
                << "Weight dimensions: " << weightDim <<nl
                << "Size dimensions: " << sizeDim
                << abort(FatalError);
        }
    }

    correct();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vdfPhaseModel::~vdfPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::vdfPhaseModel::nNodes() const
{
    return quadrature_.nodes().size();
}

Foam::tmp<Foam::volScalarField> Foam::vdfPhaseModel::volumeFraction
(
    const label nodei
) const
{
    if (nodei == -1)
    {
        tmp<volScalarField> alpha
        (
            quadrature_.moments()(volumeFractionMoment_)
        );
        if (sizeType_ == mass)
        {
            alpha.ref() /= rho();
        }
        return alpha;
    }

    if (m0VolumeFraction_)
    {
        return quadrature_.nodes()[nodei].primaryWeight();
    }

    const volVelocityNode& node = quadrature_.nodes()[nodei];
    if (momentSetType_ == numberDensityMass)
    {
        return
            node.primaryAbscissae()[sizeIndex_]
           *node.primaryWeight()
           /rho();
    }
    else if (momentSetType_ == numberDensityVolume)
    {
        return
            node.primaryAbscissae()[sizeIndex_]
           *node.primaryWeight();
    }
    else if (momentSetType_ == numberDensityLength)
    {
        return
            pow3(node.primaryAbscissae()[sizeIndex_])
           *node.primaryWeight();
    }

    NotImplemented;
    return *this;
}

Foam::tmp<Foam::volScalarField>
Foam::vdfPhaseModel::d(const label nodei) const
{
    if (sizeIndex_ == -1)
    {
        return d_;
    }

    scalar pi = Foam::constant::mathematical::pi;
    if (nodei == -1)
    {
        volScalarField m0 = quadrature_.moments()(0);
        m0.max(small);

        tmp<volScalarField> dtmp(quadrature_.moments()(sizeMoment_)/m0);

        if (sizeType_ == mass)
        {
            dtmp = Foam::max(6.0*(dtmp/rho())/pi, minD_);
        }
        else if (sizeType_ == volume)
        {
            dtmp = Foam::max(6.0*dtmp/pi, minD_);
        }
        return dtmp;
    }

    const volScalarField& size =
        quadrature_.nodes()[nodei].primaryAbscissae()[sizeIndex_];

    if (sizeType_ == length)
    {
        return Foam::max(size, minD_);
    }
    else if (sizeType_ == volume)
    {
        return Foam::max(cbrt(6.0*size/pi), minD_);
    }
    else if (sizeType_ == numberDensityLength)
    {
        return Foam::max(6.0*(size/rho())/pi, minD_);
    }

    NotImplemented;
    return size;
}

const Foam::volVectorField& Foam::vdfPhaseModel::U(const label nodei) const
{
    if (nodei == -1)
    {
        return U_;
    }
    return quadrature_.nodes()[nodei].velocityAbscissae();
}

Foam::volVectorField& Foam::vdfPhaseModel::U(const label nodei)
{
    if (nodei == -1)
    {
        return U_;
    }
    return quadrature_.nodes()[nodei].velocityAbscissae();
}

Foam::tmp<Foam::volVectorField>
Foam::vdfPhaseModel::Vs(const label nodei) const
{
    if (nodei == -1)
    {
        return tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("V", name_),
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh(),
                dimensionedVector("zero", dimVelocity, Zero)
            )
        );
    }
    return U(nodei) - U_;
}

void Foam::vdfPhaseModel::solve()
{
    populationBalance_->solve();

    const labelList& velocityIndexes =
        quadrature_.nodes()[0].velocityIndexes();

    volScalarField& alpha = *this;
    alpha = volumeFraction();

    labelList orderZero(quadrature_.momentOrders()[0].size(), 0);
    volScalarField m0(quadrature_.moments()(0));
    m0.max(residualAlpha_.value());

    forAll(velocityIndexes, cmpt)
    {
        labelList orderOne(orderZero);
        orderOne[velocityIndexes[cmpt]] = 1;

        volScalarField meanU(quadrature_.moments()(orderOne)/m0);
        U_.replace(cmpt, meanU);
    }
    phiPtr_() = fvc::flux(U_);
    alphaPhi_ = fvc::interpolate(*this)*phiPtr_();
    alphaRhoPhi_ = fvc::interpolate(rho())*alphaPhi_;

    label sizeIndex = quadrature_.nodes()[0].sizeIndex();
    if (sizeIndex != -1)
    {
        labelList orderOne(orderZero);
        orderOne[sizeIndex] = 1;
        d_ = d();
    }
}

void Foam::vdfPhaseModel::correct()
{
    quadrature_.updateMoments();

    const labelList& velocityIndexes =
        quadrature_.nodes()[0].velocityIndexes();

    labelList orderZero(quadrature_.momentOrders()[0].size(), 0);
    volScalarField m0(quadrature_.moments()(orderZero));
    m0.max(small);

    forAll(velocityIndexes, cmpt)
    {
        labelList orderOne(orderZero);
        orderOne[velocityIndexes[cmpt]] = 1;

        volScalarField meanU(quadrature_.moments()(orderOne)/m0);
        U_.replace(cmpt, meanU);

        if (computeVariance_)
        {
            forAll(velocityIndexes, cmpt2)
            {
                labelList orderOne2(orderZero);
                labelList orderTwo(orderZero);
                orderOne2[velocityIndexes[cmpt2]] = 1;
                orderTwo[velocityIndexes[cmpt]] = 1;
                orderTwo[velocityIndexes[cmpt2]] += 1;

                volScalarField meanU2(quadrature_.moments()(orderOne2)/m0);

                volScalarField coVar
                (
                    quadrature_.moments()(orderTwo)/m0 - meanU*meanU2
                );
                sigma_.ref().replace(cmpt + cmpt2*3, coVar);
            }
        }
    }

    if (computeVariance_)
    {
        Theta_.ref() = tr(sigma_())/3.0;
    }

    volScalarField& alpha = *this;
    alpha = volumeFraction();

    if (sizeIndex_ != -1)
    {
        d_ = d();
    }
}


Foam::scalar Foam::vdfPhaseModel::realizableCo() const
{
   return populationBalance_->realizableCo();
}


Foam::scalar Foam::vdfPhaseModel::CoNum() const
{
    return populationBalance_->CoNum();
}

void Foam::vdfPhaseModel::relativeTransport()
{}


void Foam::vdfPhaseModel::averageTransport()
{}


void Foam::vdfPhaseModel::solveSource()
{}

// ************************************************************************* //
