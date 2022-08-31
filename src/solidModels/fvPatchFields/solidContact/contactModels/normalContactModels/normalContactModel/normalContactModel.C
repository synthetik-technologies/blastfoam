/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

InClass
    normalContactModel

\*---------------------------------------------------------------------------*/

#include "normalContactModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(normalContactModel, 0);
defineRunTimeSelectionTable(normalContactModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

normalContactModel::normalContactModel
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const standAlonePatch& masterFaceZonePatch,
    const standAlonePatch& slaveFaceZonePatch
)
:
    name_(name),
    patch_(patch),
    masterPatchID_(masterPatchID),
    slavePatchID_(slavePatchID),
    slaveContactPointGap_
    (
        patch.boundaryMesh().mesh().boundaryMesh()[slavePatchID].nPoints(), 0.0
    )
{
    if (slavePatchID == -1)
    {
        FatalErrorIn("normalContactModel::normalContactModel")
            << "slavePatchID cannot be -1" << abort(FatalError);
    }
}


normalContactModel::normalContactModel(const normalContactModel& nm)
:
    name_(nm.name_),
    patch_(nm.patch_),
    masterPatchID_(nm.masterPatchID_),
    slavePatchID_(nm.slavePatchID_),
    slaveContactPointGap_(nm.slaveContactPointGap_)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
