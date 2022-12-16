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

#include "AitkenDisplacementRelaxation.H"
#include "valuePointPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace displacementRelaxations
{
    defineTypeNameAndDebug(Aitken, 0);
    addToRunTimeSelectionTable(displacementRelaxation, Aitken, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementRelaxations::Aitken::Aitken
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    displacementRelaxation(mesh, dict),

    relaxFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 0.01)),
    maxRelaxFactor_
    (
        dict.lookupOrDefault<scalar>("maxRelaxationFactor", 0.1)
    ),
    aitkenRelaxFactors_
    (
        coupledPatches_.size(),
        relaxFactor_
    ),

    residuals_(coupledPatches_.size()),
    prevResiduals_(coupledPatches_.size())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementRelaxations::Aitken::~Aitken()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::displacementRelaxations::Aitken::updateError
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

void Foam::displacementRelaxations::Aitken::relax
(
    const label iter,
    pointVectorField& p
)
{
    updateError(p);

    if (iter < 1)
    {
        // pointVectorField::Boundary& bp = p.boundaryFieldRef();
        // forAll(coupledPatches_, pi)
        // {
        //     const label patchi = coupledPatches_[pi];
        //     if (isA<valuePointPatchVectorField>(bp[patchi]))
        //     {
        //         valuePointPatchVectorField& pp =
        //             dynamicCast<valuePointPatchVectorField>(bp[patchi]);
        //         const valuePointPatchVectorField& ppPrev =
        //             dynamicCast<const valuePointPatchVectorField>
        //             (
        //                 p.prevIter().boundaryField()[patchi]
        //             );
        //         pp == ppPrev + relaxFactor_*residuals_[pi];
        //         pp.setInInternalField(p, pp);
        //     }
        // }
        return;
    }
    forAll(aitkenRelaxFactors_, i)
    {
        aitkenRelaxFactors_[i] =
            mag
            (
                aitkenRelaxFactors_[i]
               *(
                    gSum(prevResiduals_[i] & (residuals_[i] - prevResiduals_[i]))
                   /max
                    (
                        gSum(magSqr(residuals_[i] - prevResiduals_[i])),
                        small
                    )
                )
            );

        reduce(aitkenRelaxFactors_[i], minOp<scalar>());

        if (aitkenRelaxFactors_[i] > maxRelaxFactor_)
        {
            aitkenRelaxFactors_[i] = maxRelaxFactor_;
        }
    }

    pointVectorField::Boundary& bp = p.boundaryFieldRef();
    forAll(coupledPatches_, pi)
    {
        const label patchi = coupledPatches_[pi];
        if (isA<valuePointPatchVectorField>(bp[patchi]))
        {
            valuePointPatchVectorField& pp =
                dynamicCast<valuePointPatchVectorField>(bp[patchi]);
            const valuePointPatchVectorField& ppPrev =
                dynamicCast<const valuePointPatchVectorField>
                (
                    p.prevIter().boundaryField()[patchi]
                );
            pp == ppPrev + aitkenRelaxFactors_[pi]*residuals_[pi];
            pp.setInInternalField(p, pp);
        }
    }
}


void Foam::displacementRelaxations::Aitken::clear()
{
    displacementRelaxation::clear();
    forAll(residuals_, pi)
    {
        residuals_[pi] = Zero;
    }
}

// ************************************************************************* //
