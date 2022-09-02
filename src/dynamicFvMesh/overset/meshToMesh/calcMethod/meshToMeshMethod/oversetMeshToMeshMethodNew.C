/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "oversetMeshToMeshMethod.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::oversetMeshToMeshMethod> Foam::oversetMeshToMeshMethod::New
(
    const word& methodName,
    const polyMesh& src,
    const polyMesh& tgt
)
{
    if (debug)
    {
        Info<< "Selecting AMIMethod " << methodName << endl;
    }

    meshComponentsConstructorTable::iterator cstrIter =
        meshComponentsConstructorTablePtr_->find(methodName);

    if (cstrIter == meshComponentsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown meshToMesh type "
            << methodName << nl << nl
            << "Valid meshToMesh types are:" << nl
            << meshComponentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<oversetMeshToMeshMethod>(cstrIter()(src, tgt));
}


// ************************************************************************* //
