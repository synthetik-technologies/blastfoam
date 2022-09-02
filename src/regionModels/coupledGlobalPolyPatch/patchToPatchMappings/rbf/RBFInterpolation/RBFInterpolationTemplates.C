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

Description
    RBF interpolation templates

Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.

\*---------------------------------------------------------------------------*/

#include "RBFInterpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::RBFInterpolation::interpolate
(
    const Field<Type>& ctrlField
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>(dataPoints_.size(), pTraits<Type>::zero)
    );
    interpolate(ctrlField, tresult.ref());
    return tresult;
}


template<class Type>
void Foam::RBFInterpolation::interpolate
(
    const Field<Type>& ctrlField,
    Field<Type>& result
) const
{
    // HJ and FB (05 Jan 2009)
    // Collect the values from ALL control points to all CPUs
    // Then, each CPU will do interpolation only on local dataPoints_

    if (ctrlField.size() != controlPoints_.size())
    {
        FatalErrorInFunction
            << "Incorrect size of source field.  Size = " << ctrlField.size()
            << " nControlPoints = " << controlPoints_.size()
            << abort(FatalError);
    }


    // FB 21-12-2008
    // 1) Calculate alpha and beta coefficients using the Inverse
    // 2) Calculate displacements of internal nodes using RBF values,
    //    alpha's and beta's
    // 3) Return displacements using tresult()

    const label nControlPoints = controlPoints_.size();
    const scalarSquareMatrix& mat = this->B();

    // Determine interpolation coefficients
    Field<Type> alpha(nControlPoints, pTraits<Type>::zero);
    Field<Type> beta(4, pTraits<Type>::zero);

    for (label row = 0; row < nControlPoints; row++)
    {
        for (label col = 0; col < nControlPoints; col++)
        {
            alpha[row] += mat[row][col]*ctrlField[col];
        }
    }

    if (polynomials_)
    {
        for
        (
            label row = nControlPoints;
            row < nControlPoints + 4;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                beta[row - nControlPoints] += mat[row][col]*ctrlField[col];
            }
        }
    }

    // Evaluation
    scalar t;

    // Algorithmic improvement, Matteo Lombardi.  21/Mar/2011

    forAll (dataPoints_, flPoint)
    {
        // Cut-off function to justify neglecting outer boundary points
        t = (mag(dataPoints_[flPoint] - focalPoint_) - innerRadius_)/
            (outerRadius_ - innerRadius_);

        if (t >= 1)
        {
            // Increment is zero: w = 0
            result[flPoint] = Zero;
        }
        else
        {
            // Full calculation of weights
            scalarField weights
            (
                RBF_->weights(controlPoints_, dataPoints_[flPoint])
            );

            forAll (controlPoints_, i)
            {
                result[flPoint] += weights[i]*alpha[i];
            }

            if (polynomials_)
            {
                result[flPoint] +=
                    beta[0]
                  + beta[1]*dataPoints_[flPoint].x()
                  + beta[2]*dataPoints_[flPoint].y()
                  + beta[3]*dataPoints_[flPoint].z();
            }

            scalar w;

            if (t <= 0)
            {
                w = 1.0;
            }
            else
            {
                w = 1.0 - sqr(t)*(3.0 - 2.0*t);
            }

            result[flPoint] = w*result[flPoint];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

