/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "leastSquaresCellCellStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "SVD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(leastSquares, 0);
    addToRunTimeSelectionTable(cellCellStencil, leastSquares, mesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellCellStencils::leastSquares::stencilWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
) const
{
    // Implicit least squares weighting
    // Number of donors
    label nD = donorCcs.size();

    weights.setSize(nD);

    // List for distance vectors and LSQ weights
    List<vector> d(nD);
    scalarList LSQw(nD);

    // Sum of weights
    scalar W = 0;

    // Sum of weighted distance vectors
    vector dw(Zero);

    RectangularMatrix<scalar> A(nD, 3);

    bool shortC = false;

    // Compute distance vectors and fill rectangular matrix
    forAll(donorCcs, j)
    {
        // Neighbour weights
        d[j] = donorCcs[j] - sample;

        // Check for short-circuiting if zero distance
        // is detected with respect to any donor
        if (mag(d[j]) < ROOTVSMALL)
        {
            shortC = true;
            break;
        }

        LSQw[j] = 1.0/magSqr(d[j]);

        // T matrix
        vector wd = LSQw[j]*d[j];
        A[j][0] = wd.x();
        A[j][1] = wd.y();
        A[j][2] = wd.z();

        // Sum of weighted distance vectors
        dw += wd;

        // Sum of weights
        W += LSQw[j];
    }

    if (!shortC)
    {
        // Use Singular Value Decomposition to avoid problems
        // with 1D, 2D stencils
        SVD svd(A.T()*A, SMALL);

        // Least squares vectors
        RectangularMatrix<scalar> ATAinvAT(svd.VSinvUt()*A.T());

        scalar saveDiag(W);

        // Diagonal coefficient
        for (label i = 0; i < 3; i++)
        {
            // Get row
            scalarList Row(UList<scalar>(ATAinvAT[i], nD));

            forAll(donorCcs, j)
            {
                W -= Row[j]*LSQw[j]*dw[i];
            }
        }

        if (mag(W) < SMALL*mag(saveDiag))
        {
            shortC = true;
        }
        else
        {
            // Compute final neighbour weights with  additional scaling
            forAll(donorCcs, j)
            {
                weights[j] =
                (
                    LSQw[j]
                - ATAinvAT[0][j]*LSQw[j]*dw[0]
                - ATAinvAT[1][j]*LSQw[j]*dw[1]
                - ATAinvAT[2][j]*LSQw[j]*dw[2]
                )
            /W;
            }
        }
    }

    if (shortC)
    {
        // Matrix ill conditioned. Use straight injection from central
        // donor.
        weights = 0.0;
        weights[0] = 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::leastSquares::leastSquares
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    inverseDistance(mesh, dict, false)
{
    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::leastSquares::~leastSquares()
{}


// ************************************************************************* //
