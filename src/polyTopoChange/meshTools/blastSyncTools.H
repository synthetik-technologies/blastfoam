/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
03-10-2022  Synthetik Applied Technologies : Collected functions
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

Namespace
    Foam::syncTools

Description
    Collection of functions to do various modifications to the mesh

SourceFiles
    blastSyncTools.C
    blastSyncToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef blastSyncTools_H
#define blastSyncTools_H

#include "syncTools.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointMesh;

/*---------------------------------------------------------------------------*\
                        Namespace syncTools Declaration
\*---------------------------------------------------------------------------*/

class blastSyncTools
{

public:

     //- Template to update all processor boundaries
    template<class GeoField>
    static void correctProcessorBoundaries(const fvMesh&);

    //- Template to update all processor boundaries
    template<class Type>
    static void setInPointBoundaries(const fvMesh&);

    //- Template to update all processor boundaries
    template<class Type>
    static void correctPointBoundaries(const fvMesh&);

    //- Helper: push master point data to collocated points
    template<class Type>
    static void pushUntransformedData(const polyMesh&, Field<Type>&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blastSyncToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
