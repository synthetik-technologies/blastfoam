/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

Class
    Foam::surfacesToCell

Description
    Set a list of surfaces or a backup region

SourceFiles
    surfacesToCell.C
    surfaceToCell.C
    boxToCell.C
    cylinderToCell.C
    sphereToCell.

\*---------------------------------------------------------------------------*/

#ifndef surfacesToCell_H
#define surfacesToCell_H

#include "surfaceToCell.H"
#include "transformer.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sphericalMassToCell Declaration
\*---------------------------------------------------------------------------*/

class surfacesToCell
:
    public topoSetSource
{
// Private data

    //- Is this region a backup
    const bool isBackup_;

    //- List of trisurfaces (stored because surfaceToCell only keeps a pointer)
    PtrList<triSurface> triSurfaces_;
    PtrList<triSurfaceSearch> triSurfaceSearches_;

    //- List of surfaces
    PtrList<surfaceToCell> surfaces_;

    //- Backup topoSetSource
    autoPtr<topoSetSource> backup_;

    //- Create a transformer
    transformer createTransformer(const dictionary&) const;


public:

    //- Runtime type information
    TypeName("surfacesToCell");


    // Constructors

        //- Construct from dictionary
        surfacesToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfacesToCell();

    // Member functions

        virtual sourceType setType() const
        {
            return CELLSETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
