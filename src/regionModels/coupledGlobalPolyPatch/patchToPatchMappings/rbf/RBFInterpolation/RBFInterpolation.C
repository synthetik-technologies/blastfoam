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
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::scalarSquareMatrix& Foam::RBFInterpolation::B() const
{
    if (!BPtr_)
    {
        calcB();
    }

    return *BPtr_;
}


void Foam::RBFInterpolation::calcB() const
{
    // Determine inverse of boundary connectivity matrix
    label polySize(4);

    if (!polynomials_)
    {
        polySize = 0;
    }

    // Fill Nb x Nb matrix
    simpleMatrix<scalar> A(controlPoints_.size()+polySize);

    const label nControlPoints = controlPoints_.size();
    for (label i = 0; i < nControlPoints; i++)
    {
        scalarField weights(RBF_->weights(controlPoints_, controlPoints_[i]));

        for (label col = 0; col < nControlPoints; col++)
        {
            A[i][col] = weights[col];
        }
    }

    if (polynomials_)
    {
        for
        (
            label row = nControlPoints;
            row < nControlPoints + 1;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                A[col][row] = 1.0;
                A[row][col] = 1.0;
            }
        }

        // Fill in X components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 1;
            row < nControlPoints + 2;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                A[col][row] = controlPoints_[col].x();
                A[row][col] = controlPoints_[col].x();
            }
        }

        // Fill in Y components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 2;
            row < nControlPoints + 3;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                A[col][row] = controlPoints_[col].y();
                A[row][col] = controlPoints_[col].y();
            }
        }
        // Fill in Z components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 3;
            row < nControlPoints + 4;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                A[col][row] = controlPoints_[col].z();
                A[row][col] = controlPoints_[col].z();
            }
        }

        // Fill 4x4 zero part of matrix
        for
        (
            label row = nControlPoints;
            row < nControlPoints + 4;
            row++
        )
        {
            for
            (
                label col = nControlPoints;
                col < nControlPoints + 4;
                col++
            )
            {
                A[row][col] = 0.0;
            }
        }
    }

    // HJ and FB (05 Jan 2009)
    // Collect ALL control points from ALL CPUs
    // Create an identical inverse for all CPUs

    Info<< "Inverting RBF motion matrix" << endl;

//     BPtr_ = new scalarSquareMatrix(A.LUinvert());
    BPtr_ = new scalarSquareMatrix(SVDinv(A));
}


void Foam::RBFInterpolation::clearOut()
{
    deleteDemandDrivenData(BPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFInterpolation::RBFInterpolation
(
    const dictionary& dict,
    const vectorField& controlPoints,
    const vectorField& dataPoints
)
:
    controlPoints_(controlPoints),
    dataPoints_(dataPoints),
    RBF_(RBFFunction::New(dict)),
    BPtr_(nullptr),
    focalPoint_(dict.lookup("focalPoint")),
    innerRadius_(readScalar(dict.lookup("innerRadius"))),
    outerRadius_(readScalar(dict.lookup("outerRadius"))),
    polynomials_(dict.lookup("polynomials"))
{}


Foam::RBFInterpolation::RBFInterpolation
(
    const word& type,
    const dictionary& dict,
    const vectorField& controlPoints,
    const vectorField& dataPoints
)
:
    controlPoints_(controlPoints),
    dataPoints_(dataPoints),
    RBF_(RBFFunction::New(type, dict)),
    BPtr_(nullptr),
    focalPoint_(Zero),
    innerRadius_(0.0),
    outerRadius_(1.0),
    polynomials_(true)
{}


Foam::RBFInterpolation::RBFInterpolation
(
    const RBFInterpolation& rbf
)
:
    controlPoints_(rbf.controlPoints_),
    dataPoints_(rbf.dataPoints_),
    RBF_(rbf.RBF_->clone()),
    BPtr_(nullptr),
    focalPoint_(rbf.focalPoint_),
    innerRadius_(rbf.innerRadius_),
    outerRadius_(rbf.outerRadius_),
    polynomials_(rbf.polynomials_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBFInterpolation::~RBFInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBFInterpolation::movePoints()
{
    clearOut();
}


// ************************************************************************* //
