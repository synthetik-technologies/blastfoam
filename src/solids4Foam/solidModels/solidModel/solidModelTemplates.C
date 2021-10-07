/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "solidModel.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::solidModel::converged
(
    const int iCorr,
    const scalar solverPerfInitRes,
    const int solverPerfNIters,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const bool writeResiduals
)
{
    // We will check three residuals:
    // - relative displacement residual
    // - linear equation residual
    // - material model residual
    bool converged = false;

    // Calculate displacement residual based on the relative change of vf
    scalar denom = gMax
    (
        DimensionedField<scalar, volMesh>
        (
            mag(vf.internalField() - vf.oldTime().internalField())
        )
    );
    if (denom < SMALL)
    {
        denom = max
        (
            gMax
            (
                DimensionedField<scalar, volMesh>(mag(vf.internalField()))
            ),
            SMALL
        );
    }
    const scalar residualvf =
        gMax
        (
            DimensionedField<scalar, volMesh>
            (
                mag(vf.internalField() - vf.prevIter().internalField())
            )
        )/denom;

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at least 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol_)
    {
        if
        (
            solverPerfInitRes < solutionTol_
         && residualvf < solutionTol_
        )
        {
            if (writeResiduals)
            {
                Info<< "    Both residuals have converged" << endl;
            }
            converged = true;
        }
        else if (residualvf < alternativeTol_)
        {
            if (writeResiduals)
            {
                Info<< "    The relative residual has converged" << endl;
            }
            converged = true;
        }
        else if (solverPerfInitRes < alternativeTol_)
        {
            if (writeResiduals)
            {
                Info<< "    The solver residual has converged" << endl;
            }
            converged = true;
        }
        else
        {
            converged = false;
        }
    }

    if (!writeResiduals)
    {
        return converged;
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res, relRes, matRes, iters" << endl;
    }
    else if (iCorr % infoFrequency_ == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfInitRes
            << ", " << residualvf
            << ", " << materialResidual
            << ", " << solverPerfNIters << endl;

        if (residualFilePtr_.valid())
        {
            residualFilePtr_()
                << solverPerfInitRes << " "
                << residualvf << " "
                << materialResidual
                << endl;
        }

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr_ - 1)
    {
        maxIterReached_++;
        Warning
            << "Max iterations reached within momentum loop" << endl;
    }

    return converged;
}


// ************************************************************************* //
