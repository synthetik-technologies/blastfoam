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
    defineTypeNameAndDebug(stressCellRemovalLaw, 0);
    addToRunTimeSelectionTable
    (
        cellRemovalLaw, stressCellRemovalLaw, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stressCellRemovalLaw::stressCellRemovalLaw
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict),
    vonMises_(dict.lookupOrDefault("vonMises", true)),
    stressCrit_(dict.lookup<scalar>("criticalStress")),
    sigmaName_(dict.lookupOrDefault<word>("sigmaName", "sigma")),
    exposedPatch_(dict.lookup("exposedPatch"))
{
    if (exposedFacesPatchID() < 0)
    {
        FatalIOErrorInFunction(dict)
            << exposedPatch_ << " is not a valid patch. Valid patches are " << nl
            << mesh.boundaryMesh().names() << endl
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::stressCellRemovalLaw::~stressCellRemovalLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::stressCellRemovalLaw::cellsToRemove()
{
    // Lookup the plastic equivalent strain
    if (mesh().foundObject<volSymmTensorField>(sigmaName_))
    {
        const symmTensorField& sigma =
            mesh().lookupObject<volSymmTensorField>(sigmaName_);
        const scalarField stress
        (
            vonMises_
          ? sqrt((3.0/2.0)*magSqr(dev(sigma)))
          : mag(sigma)
        );

        labelHashSet cellsToRemove;

        forAll(stress, cellI)
        {
            if (stress[cellI] > stressCrit_)
            {
                cellsToRemove.insert(cellI);
            }
        }

        return cellsToRemove.toc();
    }
    else
    {
        return labelList(0);
    }
}


Foam::label Foam::stressCellRemovalLaw::exposedFacesPatchID()
{
    return mesh().boundaryMesh()[exposedPatch_].index();
}

// ************************************************************************* //
