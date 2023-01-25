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

#include "fixedDisplacementRelaxation.H"
#include "valuePointPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace displacementRelaxations
{
    defineTypeNameAndDebug(fixed, 0);
    addToRunTimeSelectionTable(displacementRelaxation, fixed, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementRelaxations::fixed::fixed
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    displacementRelaxation(mesh, dict),
    relaxationFactor_(coeffDict(dict).lookupOrDefault("relaxationFactor", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementRelaxations::fixed::~fixed()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::displacementRelaxations::fixed::relax
(
    const label iter,
    pointVectorField& p
)
{
    updateError(p);
    if (mesh_.relaxField(p.name()) || relaxationFactor_ < 1)
    {
        scalar f =
            mesh_.relaxField(p.name())
          ? mesh_.fieldRelaxationFactor(p.name())
          : relaxationFactor_;

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
                pp == ppPrev*(1.0 - f) + f*pp;
                pp.setInInternalField(p, pp);
            }
        }
    }

}

// ************************************************************************* //
