/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::backupSearchableSurface

Description
    Searchable surface with an optional backup

SourceFiles
    backupSearchableSurface.C

\*---------------------------------------------------------------------------*/

#ifndef backupSearchableSurface_H
#define backupSearchableSurface_H

#include "backupTopoSetSource.H"
#include "searchableSurface.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                    Class backupSearchableSurface Declaration
\*---------------------------------------------------------------------------*/

class backupSearchableSurface
{
// Private data

    //- Reference to mesh
    const polyMesh& mesh_;

    //- Stored dictionary
    dictionary dict_;

    //- Is the backup region allowed to be used
    bool allowBackup_;

    //- Searchable surface
    autoPtr<searchableSurface> surfacePtr_;

    //- Searchable surface backup
    autoPtr<searchableSurface> backupPtr_;


public:

    //- Runtime type information
    TypeName("backupSearchableSurface");

    // Constructors

        //- Copy Constructor
        backupSearchableSurface(const backupSearchableSurface&);

        //- Construct from components
        backupSearchableSurface
        (
            const word& surfaceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        // Construct from a surface
        backupSearchableSurface
        (
            searchableSurface* surface,
            const dictionary& dict
        );

    //- Destructor
    virtual ~backupSearchableSurface();

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

            autoPtr<backupSearchableSurface> operator()(Istream& is) const
            {
                word surfaceType(is);
                dictionary dict(is);
                return autoPtr<backupSearchableSurface>
                (
                    new backupSearchableSurface(surfaceType, mesh_, dict)
                );
            }
        };

    //- Clone
    autoPtr<backupSearchableSurface> clone() const
    {
        return autoPtr<backupSearchableSurface>
        (
            new backupSearchableSurface(*this)
        );
    }

    // Member Functions

        //- Return the dictionary dict
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Allow backup region to be used
        void allowBackup(const bool allow)
        {
            allowBackup_ = allow;
        }

        //- Return the searchableSurface
        const searchableSurface& surface() const
        {
            return surfacePtr_();
        }

        //- Return a subset of points that are inside
        labelList selectPoints
        (
            const pointField& points,
            const bool backupAllowed = true
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
