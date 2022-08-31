/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::searchableSurfaceToPoint

Description
    A topoSetSource to select faces based on intersection (of cell-cell
    vector) with a surface.

SourceFiles
    searchableSurfaceToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceToPoint_H
#define searchableSurfaceToPoint_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;

/*---------------------------------------------------------------------------*\
                 Class searchableSurfaceToPoint Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceToPoint
:
    public topoSetSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Surface
        autoPtr<searchableSurface> surfacePtr_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("searchableSurfaceToPoint");

    // Constructors

        //- Construct from dictionary
        searchableSurfaceToPoint
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from type name and dictionary
        searchableSurfaceToPoint
        (
            const word& surfaceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from searchable surface
        searchableSurfaceToPoint
        (
            const polyMesh& mesh,
            autoPtr<searchableSurface>& surface
        );


    //- Destructor
    virtual ~searchableSurfaceToPoint();


    // Member Functions

        //- Return the searchableSurface
        const searchableSurface& surface() const
        {
            return surfacePtr_();
        }

        virtual sourceType setType() const
        {
            return POINTSETSOURCE;
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
