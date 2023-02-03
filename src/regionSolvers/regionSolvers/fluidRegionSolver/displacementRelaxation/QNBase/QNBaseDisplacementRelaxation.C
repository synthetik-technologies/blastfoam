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

#include "QNBaseDisplacementRelaxation.H"
#include "valuePointPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace displacementRelaxations
{
    defineTypeNameAndDebug(QNBase, 0);
    addToRunTimeSelectionTable(displacementRelaxation, QNBase, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementRelaxations::QNBase::QNBase
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
    oldResiduals_(coupledPatches_.size()),

    pointsV_(coupledPatches_.size()),
    pointsW_(coupledPatches_.size()),
    times_(coupledPatches_.size())
{
    Info<< "Using " << typeName  << " relaxation with:" << nl
        << "initialRelaxationFactor: " << initRelaxFactor_ << nl
        << "nCouplingTimes: " << nOldTimes_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementRelaxations::QNBase::~QNBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::displacementRelaxations::QNBase::updateError
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


void Foam::displacementRelaxation::QNBase::updateFields
(
    const label iter,
    const pointVectorField& p
)
{
    const pointVectorField& pOld = p.oldTime();
    const pointVectorField& pPrev = p.prevIter();

    pointVectorField::Boundary& bp = p.boundaryFieldRef();
    forAll(coupledPatches_, pi)
    {
        const label patchi = coupledPatches_[pi];
        tmp<vectorField> tpp;
        tmp<vectorField> tppPrev;
        tmp<vectorField> tppOld;
        if (isA<valuePointPatchVectorField>(bp[patchi]))
        {
            tpp =
                tmp<vectorField>
                (
                    dynamicCast<valuePointPatchVectorField>(bp[patchi])
                );
            tppPrev =
                tmp<vectorField>
                (
                    dynamicCast<const valuePointPatchVectorField>
                    (
                        pPrev.boundaryField()[patchi]
                    )
                );
            tppOld =
                tmp<vectorField>
                (
                    dynamicCast<const valuePointPatchVectorField>
                    (
                        pOld.boundaryField()[patchi]
                    )
                );
        }
        else
        {
            tpp = bp[patchi].patchInternalField():
            tppPrev = pPrev.boundaryField()[patchi].patchInternalField();
            tppOld = pOld.boundaryField()[patchi].patchInternalField();
        }
        const vectorField& pp = tpp();
        const vectorField& ppPrev = tppPrev();
        const vectorField& ppOld = tppOld();

        if (iter == 1 || oldResiduals_[pi].size() != residuals_[pi].size())
        {
            oldResiduals_[pi] = residuals_[pi];
        }
        else if (iter > 1)
        {
            pointsV_[pi].append(residuals_[pi] - oldResiduals_[pi]);
            pointsW_[pi].append(pp - ppOld);
            times_[pi].append(p.time().timeIndex());
        }
    }
}


void Foam::displacementRelaxations::QNBase::clear()
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
        label startI = times.size();
        forAll(times, ti)
        {
            if (times[ti] < (mesh_.time().timeIndex() - nOldTimes_))
            {
                break;
            }
            startI = ti;
        }
        const label& oldSize = times.size();
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
    }
}

// ************************************************************************* //
