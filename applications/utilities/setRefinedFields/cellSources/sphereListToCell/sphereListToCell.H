/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::sphereListToCell

Description
    A topoSetSource to select cells based on cell centres inside sphere.

SourceFiles
    sphereListToCell.C

\*---------------------------------------------------------------------------*/

#ifndef sphereListToCell_H
#define sphereListToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sphereListToCell Declaration
\*---------------------------------------------------------------------------*/

class sphereListToCell
:
    public topoSetSource
{

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Centres
        List<vector> centres_;

        //- Radii
        List<scalar> radii_;


    // Private Member Functions

        void combine
        (
            topoSet& set,
            const vector& centre,
            const scalar& radius,
            const bool add
        ) const;


public:

    //- Runtime type information
    TypeName("sphereListToCell");


    // Constructors

        //- Construct from components
        sphereListToCell
        (
            const polyMesh& mesh,
            const List<vector>& centre,
            const List<scalar>& radius
        );

        //- Construct from dictionary
        sphereListToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        sphereListToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~sphereListToCell();


    // Member Functions

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
