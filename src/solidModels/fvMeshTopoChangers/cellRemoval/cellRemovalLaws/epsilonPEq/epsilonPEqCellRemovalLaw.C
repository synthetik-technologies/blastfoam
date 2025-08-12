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

#include "epsilonPEqCellRemovalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "removeCells.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellRemovalLaws
{
    defineTypeNameAndDebug(epsilonPEq, 0);
    addToRunTimeSelectionTable
    (
        cellRemovalLaw, epsilonPEq, dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellRemovalLaws::epsilonPEq::epsilonPEq
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict),
    epsilonPEqCrit_(readScalar(dict.lookup("epsilonPEqCritical"))),
    epsilonPEqName_(dict.lookupOrDefault<word>("epsilonPEqName", "epsilonPEq"))
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::cellRemovalLaws::epsilonPEq::~epsilonPEq()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelHashSet Foam::cellRemovalLaws::epsilonPEq::cellsToRemove()
{
    // Lookup the plastic equivalent strain
    const volScalarField& epsilonPEq =
        mesh_.lookupObject<volScalarField>(epsilonPEqName_);

    labelHashSet cellsToRemoveSet;
    forAll(epsilonPEq, cellI)
    {
        if (epsilonPEq[cellI] > epsilonPEqCrit_)
        {
            cellsToRemoveSet.insert(cellI);
        }
    }

    return cellsToRemoveSet;
}


// ************************************************************************* //
