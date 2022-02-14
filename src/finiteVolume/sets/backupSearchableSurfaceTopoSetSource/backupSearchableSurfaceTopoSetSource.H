/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::backupSearchableSurfaceTopoSetSource

Description
    Runtime selectable topoSetSource with a backup region

SourceFiles
    backupSearchableSurfaceTopoSetSource.C

\*---------------------------------------------------------------------------*/

#ifndef backupSearchableSurfaceTopoSetSource_H
#define backupSearchableSurfaceTopoSetSource_H

#include "backupTopoSetSource.H"
#include "searchableSurface.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                    Class backupSearchableSurfaceTopoSetSource Declaration
\*---------------------------------------------------------------------------*/

class backupSearchableSurfaceTopoSetSource
:
    public backupTopoSetSource
{
// Private data

    //- Searchable surface
    UautoPtr<const searchableSurface> surfacePtr_;


public:

    //- Runtime type information
    TypeName("backupSearchableSurfaceTopoSetSource");

    // Constructors

        //- Construct from components
        backupSearchableSurfaceTopoSetSource
        (
            const word& topoSetSourceType,
            const polyMesh& mesh,
            const dictionary& dict
        );


        //- Class used for the read-construction of
        //  PtrLists of backupTopoSetSource
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<backupSearchableSurfaceTopoSetSource> operator()
            (
                Istream& is
            ) const
            {
                word topoSetSourceType(is);
                dictionary dict(is);
                return autoPtr<backupSearchableSurfaceTopoSetSource>
                (
                    new backupSearchableSurfaceTopoSetSource
                    (
                        topoSetSourceType,
                        mesh_,
                        dict
                    )
                );
            }
        };

    //- Destructor
    virtual ~backupSearchableSurfaceTopoSetSource();



    //- Clone
    autoPtr<backupSearchableSurfaceTopoSetSource> clone() const
    {
        return autoPtr<backupSearchableSurfaceTopoSetSource>
        (
            new backupSearchableSurfaceTopoSetSource(*this)
        );
    }

    // Member Functions

        //- Return the searchableSurface
        const searchableSurface& surface() const
        {
            return surfacePtr_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
