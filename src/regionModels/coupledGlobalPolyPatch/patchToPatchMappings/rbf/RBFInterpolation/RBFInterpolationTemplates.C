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
    const Field<Type>& fromField
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>(dataPoints_.size(), pTraits<Type>::zero)
    );
    interpolate(fromField, tresult.ref());
    return tresult;
}


// template<class Type>
// Foam::tmp<Foam::Field<Type>> Foam::RBFInterpolation::interpolate2
// (
//     const Field<Type>& fromField
// ) const
// {
//     tmp<Field<Type> > tresult
//     (
//         new Field<Type>(dataPoints_.size(), pTraits<Type>::zero)
//     );
//     interpolate2(fromField, tresult.ref());
//     return tresult;
// }


/*
* Compute interpolation matrix and directly interpolate the values.
* The algorithms solves for the coefficients, and explicitly
* uses the coefficients to interpolate the data to the new positions.
*/
template<class Type>
void Foam::RBFInterpolation::interpolate
(
    const Field<Type> & fromField,
    Field<Type> & toField
) const
{
    const CVectorMatrix<Type> values
    (
        fromField.cdata()->v_,
        fromField.size(),
        pTraits<Type>::nComponents
    );
    VectorMatrix<Type> valuesInterpolation
    (
        toField.begin()->v_,
        toField.size(),
        pTraits<Type>::nComponents
    );

    for (label cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
    {
        Eigen::VectorXd valuesI(fromField.size() + polySize_);
        valuesI.setZero();
        forAll(fromField, i)
        {
            valuesI(i) = values(i, cmpti);
        }

        if (consistent_)
        {
            auto valuesInterpolationI = valuesInterpolation.col(cmpti);

            Eigen::VectorXd polynomialContribution;
            // Solve polynomial QR and subtract it from the input data
            if (separate_)
            {
                polynomialContribution = QR().solve(valuesI);
                valuesI -= (Q()*polynomialContribution);
            }

            // Integrated polynomial (and separated)
            Eigen::VectorXd p;
            if (RBF_->positiveDefinite())
            {
                p = decompLLT().solve(valuesI);
            }
            else
            {
                p = decompQR().solve(valuesI);
            }
            valuesInterpolationI = A()*p;

            // Add the polynomial part again for separated polynomial
            if (separate_)
            {
                valuesInterpolationI += (V()*polynomialContribution);
            }
        }
        else
        {
            Eigen::VectorXd Au = A().transpose()*valuesI;
            Eigen::VectorXd valuesInterpolationI;
            if (RBF_->positiveDefinite())
            {
                valuesInterpolationI = decompLLT().solve(Au);
            }
            else
            {
                valuesInterpolationI = decompQR().solve(Au);
            }

            if (separate_)
            {
                Eigen::VectorXd epsilon = V().transpose()*valuesI;

                // epsilon = Q^T * mu - epsilon (tau in the PETSc impl)
                epsilon -= Q().transpose()*valuesInterpolationI;

                // out  = out - solveTranspose tau (sigma in the PETSc impl)
                // Newer version of eigen provide the solve() for transpose()
                // matrix decopmositions
#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
                valuesInterpolationI -=
                    static_cast<Eigen::VectorXd>(QR().transpose().solve(-epsilon));
#else
                // Backwards compatible version
                Eigen::VectorXd sigma(Q().rows());
                const Eigen::Index nonzeroPivots = QR().nonzeroPivots();

                if (nonzeroPivots == 0)
                {
                    sigma.setZero();
                }
                else
                {
                    Eigen::VectorXd c
                    (
                        QR().colsPermutation().transpose()*(-epsilon)
                    );

                    QR().matrixQR().topLeftCorner
                    (
                        nonzeroPivots,
                        nonzeroPivots
                    ).template triangularView<Eigen::Upper>().transpose().conjugate().solveInPlace
                    (
                        c.topRows(nonzeroPivots)
                    );

                    sigma.topRows(nonzeroPivots) = c.topRows(nonzeroPivots);
                    sigma.bottomRows(QR().rows() - nonzeroPivots).setZero();

                    sigma.applyOnTheLeft(QR().householderQ().setLength(nonzeroPivots));
                    valuesInterpolationI -= sigma;
                }
#endif
                forAll(toField, i)
                {
                    valuesInterpolation(i, cmpti) = valuesInterpolationI(i);
                }
            }
        }
    }
}


// /*
// * This function is only called by the RBFCoarsening class.
// * It is assumed that the polynomial term is included in the
// * interpolation, and that the fullPivLu decomposition is
// * used to solve for the coefficients B.
// */
// template<class Type>
// void Foam::RBFInterpolation::interpolate2
// (
//     const Field<Type>& fromField,
//     Field<Type>& toField
// ) const
// {
//     const CVectorMatrix<Type> values
//     (
//         fromField.cdata()->v_,
//         fromField.size(),
//         pTraits<Type>::nComponents
//     );
//     VectorMatrix<Type> valuesInterpolation
//     (
//         toField.begin()->v_,
//         toField.size(),
//         pTraits<Type>::nComponents
//     );
//
//     for (label cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
//     {
//         auto valuesI = values.col(cmpti);
//         auto valuesInterpolationI = valuesInterpolation.col(cmpti);
//
//         if (polynomials_)
//         {
//             Matrix valuesLU
//             (
//                 values.rows() + values.cols() + 1,
//                 values.cols()
//             );
//             valuesLU.setZero(); // initialize all values zero
//             valuesLU.topLeftCorner(values.rows(), values.cols()) = valuesI;
//             valuesInterpolationI = Phi()*lu().solve(valuesLU);
//         }
//         else
//         {
//             valuesInterpolationI = Phi()*lu().solve(valuesI);
//         }
//     }
// }



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

