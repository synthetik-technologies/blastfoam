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

#include "IterativeSolidModel.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidModel>
template<class Type>
bool Foam::IterativeSolidModel<SolidModel>::converged
(
    const label iCorr,
    const scalar solverPerfInitRes,
    const label solverPerfNIters,
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
    scalar refResidual = great;

    // If this is the delta displacement field use the max of the current
    // and old time as the reference
    if (vf.name() == this->DD().name())
    {
        refResidual =
            max
            (
                gMax(mag(vf.primitiveField())),
                gMax(mag(vf.oldTime().primitiveField()))
            );
    }
    else
    {
        refResidual =
            gMax
            (
                mag(vf.primitiveField() - vf.oldTime().primitiveField())
            );
    }
    scalar residual =
        gMax(mag(vf.primitiveField() - vf.prevIter().primitiveField()));
    scalar relResidual = 0.0;
    if (refResidual > SMALL)
    {
        relResidual = residual/refResidual;
    }

    // Calculate material residual
    const scalar materialResidual = this->mechanical().residual();
    const scalar materialRelResidual = this->mechanical().relResidual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at least 1 outer iteration and the material law must be converged
    if
    (
        iCorr > 1
     && (
            materialResidual < materialTol_
         || materialRelResidual < materialRelTol_
        )
    )
    {
        if
        (
            solverPerfInitRes < solutionTol_
         && (relResidual < relTol_ || residual < tolerance_)
        )
        {
            if (writeResiduals)
            {
                Info<< nl << "    Both residuals have converged" << endl;
            }
            converged = true;
        }
        else if (residual < tolerance_ || relResidual < alternativeTol_)
        {
            if (writeResiduals)
            {
                Info<< nl
                    << "    The residual has converged to the alternative tolerance"
                    << endl;
            }
            converged = true;
        }
        else if (solverPerfInitRes < alternativeTol_)
        {
            if (writeResiduals)
            {
                Info<< nl
                    << "    The solver residual has converged to the alternative "
                    << "tolerance" << endl;
            }
            converged = true;
        }
        else
        {
            converged = false;
        }
    }
    reduce(converged, andOp<bool>());

    if (!writeResiduals)
    {
        return converged;
    }

    // Print residual information
    if (iCorr == 0)
    {
        int width = Info().precision() + 6;
        Info<< "    "
            << setf(ios_base::left)
            << setw(10) << "Corr"
            << setw(20) << "solutionResidual"
            << setw(2*width + 3) << "residual (abs/rel)"
            << setw(2*width + 3) << "materialResidual (abs/rel)"
            << setw(5) << "iters"
            << endl;
    }

    if (iCorr % infoFrequency_ == 0 || converged)
    {
        int width = Info().precision() + 6;
        Info<< "    "
            << setf(ios_base::left)
            << setw(10) << iCorr
            << setw(20) << solverPerfInitRes
            << setw(2*width + 3)
                << word(Foam::name(residual) + " / " + Foam::name(relResidual))
            << setw(2*width + 3)
                << word(Foam::name(materialResidual) + " / " + Foam::name(materialRelResidual))
            << setw(5) << solverPerfNIters << endl;

        if (residualFilePtr_.valid())
        {
            residualFilePtr_()
                << solverPerfInitRes << token::SPACE
                << residual << token::SPACE
                << relResidual << token::SPACE
                << materialResidual << token::SPACE
                << materialRelResidual << token::SPACE
                << solverPerfNIters
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
        Warning<< "Max iterations reached within momentum loop" << endl;
    }

    return converged;
}


// ************************************************************************* //
