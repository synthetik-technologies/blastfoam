/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.

\*---------------------------------------------------------------------------*/

#include "RBFInterpolation.H"
#include "QRMatrix.H"
#include "demandDrivenData.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBFInterpolation, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::RBFInterpolation::fillPolynomialEntries
(
    Matrix& M,
    const label start,
    const UList<vector>& pts
) const
{
    for (label i = 0; i < pts.size(); i++)
    {
        M(i, start) = 1.0;
    }

    label k = 0;
    for (label cmpti = 0; cmpti < 3; cmpti++)
    {
        if (validDirs_[cmpti])
        {
            for (label i = 0; i < pts.size(); i++)
            {
                M(i, start + 1 + k) = pts[i][cmpti];
            }
            k++;
        }
    }
}


void Foam::RBFInterpolation::calc() const
{
    DebugInfo<< "Computing decomposition matricies" <<endl;

    // Sizes
    const label n = controlPoints_.size() + polySize_;

    // Evaluate radial basis functions for matrix H
    Matrix CLU(n, n);
    CLU.setZero();

    const scalar rCutOffSqr = sqr(RBF_->cutOffRadius());

    // RBF function evaluation
    forAll(controlPoints_, i)
    {
        for (label j = i; j < controlPoints_.size(); j++)
        {
            scalar rSqr = magSqr(controlPoints_[i] - controlPoints_[j]);
            CLU(i, j) = rSqr < rCutOffSqr ? RBF_->evaluate(sqrt(rSqr)) : 0;
        }
    }

    // Include polynomial contributions
    if (polynomials_)
    {
        fillPolynomialEntries
        (
            CLU,
            controlPoints_.size(),
            controlPoints_
        );
    }
    CLU.triangularView<Eigen::Lower>() = CLU.transpose();

    if (decompLLTPtr_ || decompQRPtr_)
    {
        FatalErrorInFunction
            << "Decomposition matrix already set" << endl
            << abort(FatalError);
    }

    bool sucessful = false;
    if (RBF_->positiveDefinite())
    {
        decompLLTPtr_ = new LLTMatrix(CLU.llt());
        sucessful =
            decompLLTPtr_->info() == Eigen::ComputationInfo::Success;
    }
    else
    {
        decompQRPtr_ = new QRMatrix(CLU.colPivHouseholderQr());
        sucessful = decompQRPtr_->isInvertible();
    }
    if (!sucessful)
    {
        FatalErrorInFunction
             << "Interpolation matrix is not invertable" << endl
             << abort(FatalError);
    }

    if (separate_)
    {
        if (QPtr_ || VPtr_ || QRPtr_)
        {
            FatalErrorInFunction
                << "Separate matrices already set" << endl
                << abort(FatalError);
        }
        QPtr_ = new Matrix(controlPoints_.size(), 4 - nDead_);
        fillPolynomialEntries(*QPtr_, 0, controlPoints_);

        VPtr_ = new Matrix(dataPoints_.size(), 4 - nDead_);
        fillPolynomialEntries(*VPtr_, 0, dataPoints_);

        QRPtr_ = new QRMatrix(QPtr_->colPivHouseholderQr());
    }

    DebugInfo<< "Done computing decomposition matrices" << endl;
}


void Foam::RBFInterpolation::calcA() const
{
    if (APtr_)
    {
        FatalErrorInFunction
            << "A matrix already set" << endl
            << abort(FatalError);
    }
    DebugInfo<< "Computing A matrix" << endl;

    APtr_ =
        new Matrix
        (
            dataPoints_.size(),
            controlPoints_.size() + polySize_
        );
    Matrix& A = *APtr_;
    A.setZero();

    const scalar cutOff = RBF_->cutOffRadius();

    // Evaluate A which contains the evaluation of the radial basis function
    forAll(dataPoints_, i)
    {
        forAll(controlPoints_, j)
        {
            scalar r = mag(controlPoints_[j] - dataPoints_[i]);
            if (r < cutOff)
            {
                A(i, j) = RBF_->evaluate(r);
            }
        }
    }

    // Include polynomial contributions in matrix A
    if (polynomials_)
    {
        fillPolynomialEntries
        (
            A,
            controlPoints_.size(),
            dataPoints_
        );
    }
    DebugInfo<< "Finished computing A matrix" << endl;
}


void Foam::RBFInterpolation::calcHhat() const
{
    DebugInfo<< "Computing conservative interpolation matrix"<<endl;

    // Sizes
    const label n = controlPoints_.size() + polySize_;
    Matrix I(Matrix::Identity(n, n));

    if (RBF_->positiveDefinite())
    {
        HhatPtr_ = new Matrix(A()*decompLLT().solve(I));
    }
    else
    {
        HhatPtr_ = new Matrix(A()*decompQR().solve(I));
    }

    DebugInfo<< "Finished computing conservative interpolation matrix" << endl;
}


const Foam::RBFInterpolation::Matrix& Foam::RBFInterpolation::A() const
{
    if (!APtr_)
    {
        calcA();
    }

    return *APtr_;
}


const Foam::RBFInterpolation::LLTMatrix&
Foam::RBFInterpolation::decompLLT() const
{
    if (!decompLLTPtr_)
    {
        calc();
    }

    return *decompLLTPtr_;
}


const Foam::RBFInterpolation::QRMatrix&
Foam::RBFInterpolation::decompQR() const
{
    if (!decompQRPtr_)
    {
        calc();
    }

    return *decompQRPtr_;
}


const Foam::RBFInterpolation::Matrix& Foam::RBFInterpolation::Hhat() const
{
    if (!HhatPtr_)
    {
        calcHhat();
    }

    return *HhatPtr_;
}


const Foam::RBFInterpolation::Matrix& Foam::RBFInterpolation::Q() const
{
    if (!QPtr_)
    {
        calc();
    }

    return *QPtr_;
}


const Foam::RBFInterpolation::Matrix& Foam::RBFInterpolation::V() const
{
    if (!VPtr_)
    {
        calc();
    }

    return *VPtr_;
}


const Foam::RBFInterpolation::QRMatrix& Foam::RBFInterpolation::QR() const
{
    if (!QRPtr_)
    {
        calc();
    }

    return *QRPtr_;
}

void Foam::RBFInterpolation::clearOut()
{
    deleteDemandDrivenData(APtr_);
    deleteDemandDrivenData(decompLLTPtr_);
    deleteDemandDrivenData(decompQRPtr_);
    deleteDemandDrivenData(HhatPtr_);
    deleteDemandDrivenData(QPtr_);
    deleteDemandDrivenData(VPtr_);
    deleteDemandDrivenData(QRPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFInterpolation::RBFInterpolation
(
    const dictionary& dict,
    const vectorField& controlPoints,
    const vectorField& dataPoints,
    const patchToPatchMapping::Mapping mapping
)
:
    mapping_(mapping),
    controlPoints_(controlPoints),
    dataPoints_(dataPoints),
    RBF_(RBFFunction::New(dict)),

    positions_
    (
        const_cast<scalar*>(controlPoints_.begin()->v_),
        controlPoints_.size(),
        vector::nComponents
    ),
    positionsInterpolation_
    (
        const_cast<scalar*>(dataPoints_.begin()->v_),
        dataPoints_.size(),
        vector::nComponents
    ),
    APtr_(nullptr),
    decompLLTPtr_(nullptr),
    decompQRPtr_(nullptr),
    HhatPtr_(nullptr),
    QPtr_(nullptr),
    VPtr_(nullptr),
    QRPtr_(nullptr),
    separate_(dict.lookup<bool>("separate")),
    polynomials_
    (
        RBF_->positiveDefinite() || separate_
      ? false
      : dict.lookup<bool>("polynomials")
    ),
    validDirs_(1, 1, 1),
    nDead_(0),
    polySize_(0)
{
    DebugInfo<< "Using "
        << patchToPatchMapping::MappingTypeNames[mapping_]
        << " mapping " << endl;
    if (polynomials_)
    {
        DebugInfo<< "Using polynomial terms" << endl;
        polySize_ = 4;
        boundBox bb(controlPoints_);
        for (label cmpti = 0; cmpti < 3; cmpti++)
        {
            if (mag(bb.min()[cmpti] - bb.max()[cmpti]) < small)
            {
                validDirs_[cmpti] = 0;
                nDead_++;
                polySize_--;
            }
        }
    }
}


Foam::RBFInterpolation::RBFInterpolation
(
    const word& type,
    const dictionary& dict,
    const vectorField& controlPoints,
    const vectorField& dataPoints,
    const patchToPatchMapping::Mapping mapping
)
:
    mapping_(mapping),
    controlPoints_(controlPoints),
    dataPoints_(dataPoints),
    RBF_(RBFFunction::New(type, dict)),

    positions_
    (
        const_cast<scalar*>(controlPoints_.begin()->v_),
        controlPoints_.size(),
        vector::nComponents
    ),
    positionsInterpolation_
    (
        const_cast<scalar*>(dataPoints_.begin()->v_),
        dataPoints_.size(),
        vector::nComponents
    ),
    APtr_(nullptr),
    decompLLTPtr_(nullptr),
    decompQRPtr_(nullptr),
    HhatPtr_(nullptr),
    QPtr_(nullptr),
    VPtr_(nullptr),
    QRPtr_(nullptr),
    separate_(dict.lookupOrDefault("separate", false)),
    polynomials_
    (
        RBF_->positiveDefinite() || separate_
      ? false
      : dict.lookupOrDefault("polynomials", true)
    ),
    validDirs_(1, 1, 1),
    nDead_(0),
    polySize_(0)
{
     DebugInfo<< "Using "
        << patchToPatchMapping::MappingTypeNames[mapping_]
        << " mapping " << endl;
    if (polynomials_)
    {
        DebugInfo<< "Using polynomial terms" << endl;
        polySize_ = 4;
        boundBox bb(controlPoints_);
        for (label cmpti = 0; cmpti < 3; cmpti++)
        {
            if (mag(bb.min()[cmpti] - bb.max()[cmpti]) < small)
            {
                validDirs_[cmpti] = 0;
                nDead_++;
                polySize_--;
            }
        }
    }
}


Foam::RBFInterpolation::RBFInterpolation
(
    const RBFInterpolation& rbf
)
:
    mapping_(rbf.mapping_),
    controlPoints_(rbf.controlPoints_),
    dataPoints_(rbf.dataPoints_),
    RBF_(rbf.RBF_->clone()),

    positions_
    (
        const_cast<scalar*>(controlPoints_.begin()->v_),
        controlPoints_.size(),
        vector::nComponents
    ),
    positionsInterpolation_
    (
        const_cast<scalar*>(dataPoints_.begin()->v_),
        dataPoints_.size(),
        vector::nComponents
    ),
    APtr_(nullptr),
    decompLLTPtr_(nullptr),
    decompQRPtr_(nullptr),
    HhatPtr_(nullptr),
    QPtr_(nullptr),
    VPtr_(nullptr),
    QRPtr_(nullptr),
    separate_(rbf.separate_),
    polynomials_(rbf.polynomials_),
    validDirs_(1, 1, 1),
    nDead_(0),
    polySize_(0)
{
     DebugInfo<< "Using "
        << patchToPatchMapping::MappingTypeNames[mapping_]
        << " mapping " << endl;
    if (polynomials_)
    {
        DebugInfo<< "Using polynomial terms" << endl;
        polySize_ = 4;
        boundBox bb(controlPoints_);
        for (label cmpti = 0; cmpti < 3; cmpti++)
        {
            if (mag(bb.min()[cmpti] - bb.max()[cmpti]) < small)
            {
                validDirs_[cmpti] = 0;
                nDead_++;
                polySize_--;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBFInterpolation::~RBFInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
const Foam::scalar* Foam::RBFInterpolation::getData(const UList<scalar>& lst)
{
    return lst.cdata();
}

template<>
Foam::scalar* Foam::RBFInterpolation::getData(UList<scalar>& lst)
{
    return lst.data();
}


void Foam::RBFInterpolation::movePoints()
{
    clearOut();
}


// ************************************************************************* //
