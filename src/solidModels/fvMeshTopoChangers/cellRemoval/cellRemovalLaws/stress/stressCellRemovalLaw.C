/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "stressCellRemovalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "removeCells.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellRemovalLaws
{
    defineTypeNameAndDebug(stress, 0);
    addToRunTimeSelectionTable
    (
        cellRemovalLaw, stress, dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellRemovalLaws::stress::stress
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict),
    vonMises_(dict.lookupOrDefault("vonMises", true)),
    stressCrit_(dict.lookup<scalar>("criticalStress")),
    sigmaName_(dict.lookupOrDefault<word>("sigmaName", "sigma"))
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::cellRemovalLaws::stress::~stress()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelHashSet Foam::cellRemovalLaws::stress::cellsToRemove()
{
    // Lookup the plastic equivalent strain
    const symmTensorField& sigma =
        mesh_.lookupObject<volSymmTensorField>(sigmaName_);
    const scalarField stress
    (
        vonMises_
      ? sqrt((3.0/2.0)*magSqr(dev(sigma)))
      : mag(sigma)
    );

    labelHashSet cellsToRemoveSet;
    forAll(stress, cellI)
    {
        if (stress[cellI] > stressCrit_)
        {
            cellsToRemoveSet.insert(cellI);
        }
    }

    return cellsToRemoveSet;
}

// ************************************************************************* //
