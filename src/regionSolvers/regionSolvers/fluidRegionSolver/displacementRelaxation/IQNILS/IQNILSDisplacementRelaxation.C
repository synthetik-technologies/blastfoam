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

#include "IQNILSDisplacementRelaxation.H"
#include "valuePointPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace displacementRelaxations
{
    defineTypeNameAndDebug(IQNILS, 0);
    addToRunTimeSelectionTable(displacementRelaxation, IQNILS, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementRelaxations::IQNILS::IQNILS
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    displacementRelaxation(mesh, dict),

    initRelaxFactor_
    (
        coeffDict(dict).lookup<scalar>("initialRelaxationFactor")
    ),
    nOldTimes_(coeffDict(dict).lookup<label>("nCouplingTimes")),

    residuals_(coupledPatches_.size()),
    prevResiduals_(coupledPatches_.size()),

    refResiduals_(coupledPatches_.size()),
    refDisplacements_(coupledPatches_.size()),

    pointsV_(coupledPatches_.size()),
    pointsW_(coupledPatches_.size()),
    times_(coupledPatches_.size())
{
    Info<< "Using " << typeName  << " relaxation with:" << nl
        << "initialRelaxationFactor: " << initRelaxFactor_ << nl
        << "nCouplingTimes: " << nOldTimes_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementRelaxations::IQNILS::~IQNILS()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::displacementRelaxations::IQNILS::updateError
(
    const pointVectorField& p
)
{
    error_ = Zero;
    forAll(coupledPatches_, pi)
    {
        const label patchi = coupledPatches_[pi];
        const pointPatchVectorField& pp = p.boundaryField()[patchi];
        const pointPatchVectorField& ppPrev = p.prevIter().boundaryField()[patchi];
        if (isA<valuePointPatchVectorField>(pp))
        {
            prevResiduals_[pi].transfer(residuals_[pi]);
            residuals_[pi] =
                dynamicCast<const valuePointPatchVectorField>(pp)
              - dynamicCast<const valuePointPatchVectorField>(ppPrev);
            error_ += sum(magSqr(residuals_[pi]));
        }
    }
    reduce(error_, sumOp<scalar>());
    error_ = sqrt(error_);

    if (initialError_ < 0)
    {
        initialError_ = error_;
    }
}

void Foam::displacementRelaxations::IQNILS::relax
(
    const label iter,
    pointVectorField& p
)
{
    updateError(p);

    const pointVectorField& pPrev = p.prevIter();
    pointVectorField::Boundary& bp = p.boundaryFieldRef();

    forAll(coupledPatches_, pi)
    {
        const label patchi = coupledPatches_[pi];
        if (!isA<valuePointPatchVectorField>(bp[patchi]))
        {
            continue;
        }
        valuePointPatchVectorField& pp =
            dynamicCast<valuePointPatchVectorField>(bp[patchi]);
        const valuePointPatchVectorField& ppPrev =
            dynamicCast<const valuePointPatchVectorField>
            (
                pPrev.boundaryField()[patchi]
            );

        if (iter == 0)
        {
            refResiduals_[pi] = residuals_[pi];
            refDisplacements_[pi] = pp;

            pp == ppPrev + initRelaxFactor_*residuals_[pi];
            pp.setInInternalField(p, pp);
        }
        else
        {
            pointsV_[pi].append(residuals_[pi] - refResiduals_[pi]);
            pointsW_[pi].append(pp - refDisplacements_[pi]);
            times_[pi].append(p.time().timeIndex());

            label n = pointsV_[pi].size();
            scalarSquareMatrix R(n, 0.0);
            scalarField C(n, 0.0);
            scalarField RColSum(n, 0.0);
            List<vectorField> Q(n);

            for (label i = 0; i < n; i++)
            {
                Q[i] = pointsV_[pi][n-1-i];
            }

            for (label i = 0; i < n; i++)
            {
                R[i][i] = sqrt(gSum(Q[i] & Q[i]));
                Q[i] /= max(R[i][i], small);

                for (label j = i+1; j < n; j++)
                {
                    R[i][j] = gSum(Q[i] & Q[j]);
                    Q[j] -= R[i][j]*Q[i];
                }

                C[i] = gSum(Q[i] & -residuals_[pi]);
            }

            for (label j = 0; j < n; j++)
            {
                RColSum[j] = 0.0;
                for (label i = 0; i < j+1; i++)
                {
                    RColSum[j] = mag(R[i][j]);
                }
            }

            scalar epsilon = 1e-10*max(RColSum);

            for (label i = 0; i < n; i++)
            {
                if (mag(R[i][i]) > epsilon)
                {
                    for (label j = i+1; j < n; j++)
                    {
                        R[i][j] /= R[i][i];
                    }
                    C[i] /= R[i][i];
                    R[i][i] = 1.0;
                }
            }

            for (label j = n-1; j >= 0; j--)
            {
                if (mag(R[j][j]) > epsilon)
                {
                    for (label i = 0; i < j; i++)
                    {
                        C[i] -= C[j]*R[i][j];
                    }
                }
                else
                {
                    C[j] = 0.0;
                }
            }
            vectorField newDisp(pp);
            forAll(pointsW_[pi], i)
            {
                newDisp += pointsW_[pi][i]*C[n-1-i];
            }
            Info<<gMaxMagSqr(pp)<<" "<<gMaxMagSqr(newDisp)<<endl;
            pp == newDisp;
            pp.setInInternalField(p, pp);
        }
    }
}


void Foam::displacementRelaxations::IQNILS::clear()
{
    displacementRelaxation::clear();
    forAll(residuals_, pi)
    {
        residuals_[pi] = Zero;
    }

    if (!nOldTimes_)
    {
        forAll(times_, pi)
        {
            pointsV_[pi].clear();
            pointsW_[pi].clear();
            times_[pi].clear();
        }
        return;
    }

    forAll(times_, pi)
    {
        DynamicList<scalar>& times = times_[pi];
        Info<<"times: "<<times<<endl;
        label startI = times.size();

        forAll(times, ti)
        {
            Info<<(mesh_.time().timeIndex() - times[ti])<<" "<<nOldTimes_<<endl;
            if ((mesh_.time().timeIndex() - times[ti]) < nOldTimes_)
            {
                startI = ti;
                break;
            }
        }

        const label oldSize = times.size();
        if (startI > 0)
        {
            for (label ti = 0; ti < oldSize-startI; ti++)
            {
                pointsV_[pi][ti].transfer(pointsV_[pi][ti+startI]);
                pointsW_[pi][ti].transfer(pointsW_[pi][ti+startI]);
                times[ti] = times[ti+startI];
            }
            pointsV_[pi].setSize(oldSize-startI);
            pointsW_[pi].setSize(oldSize-startI);
            times.setSize(oldSize-startI);
        }
        Info<<times<<endl<<endl;
    }
}

// ************************************************************************* //
