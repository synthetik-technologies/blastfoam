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

#include "positionCellRemovalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "cellSet.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(positionCellRemovalLaw, 0);
    addToRunTimeSelectionTable
    (
        cellRemovalLaw, positionCellRemovalLaw, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::positionCellRemovalLaw::positionCellRemovalLaw
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict),
    exposedPatch_(dict.lookup("exposedPatch"))
{
    if (exposedFacesPatchID() < 0)
    {
        FatalIOErrorInFunction(dict)
            << exposedPatch_ << " is not a valid patch. Valid patches are " << nl
            << mesh.boundaryMesh().names() << endl
            << abort(FatalIOError);
    }

    PtrList<entry> sourceDicts(dict.lookup("cellsToRemove"));
    sources_.setSize(sourceDicts.size());
    times_.setSize(sourceDicts.size());

    forAll(sourceDicts, i)
    {
        sources_.set
        (
            i,
            topoSetSource::New
            (
                sourceDicts[i].keyword(),
                mesh,
                sourceDicts[i].dict()
            )
        );
        times_[i] = sourceDicts[i].dict().lookup<scalar>("time");
    }

    labelList order;
    sortedOrder(times_, order);
    inplaceReverseList(order);

    times_ = scalarList(times_, order);
    sources_.shuffle(order);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::positionCellRemovalLaw::~positionCellRemovalLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::positionCellRemovalLaw::cellsToRemove()
{
    if (times_.size() && mesh().time().value() > times_.last())
    {
        cellSet cset(mesh(), "cset", mesh().nCells()/10+1);
        sources_.last().applyToSet(topoSetSource::NEW, cset);

        times_.setSize(times_.size()-1);
        sources_.setSize(sources_.size()-1);
        return cset.toc();
    }
    return labelList::null();

}


Foam::label Foam::positionCellRemovalLaw::exposedFacesPatchID()
{
    return mesh().boundaryMesh()[exposedPatch_].index();
}

// ************************************************************************* //
