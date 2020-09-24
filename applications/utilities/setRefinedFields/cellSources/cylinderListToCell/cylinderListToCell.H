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
    Foam::cylinderListToCell

Description
    A topoSetSource to select cells based on cell centres inside a list
    of cylinders.

SourceFiles
    cylinderListToCell.C

\*---------------------------------------------------------------------------*/

#ifndef cylinderListToCell_H
#define cylinderListToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cylinderListToCell Declaration
\*---------------------------------------------------------------------------*/

class cylinderListToCell
:
    public topoSetSource
{

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- First points on cylinder axis
        List<vector> p1s_;

        //- Second points on cylinder axis
        List<vector> p2s_;

        //- Radii
        List<scalar> radii_;


    // Private Member Functions

        void combine
        (
            topoSet& set,
            const vector& p1,
            const vector& p2,
            const scalar& radius,
            const bool add
        ) const;


public:

    //- Runtime type information
    TypeName("cylinderListToCell");


    // Constructors

        //- Construct from components
        cylinderListToCell
        (
            const polyMesh& mesh,
            const List<vector>& p1s,
            const List<vector>& p2s,
            const List<scalar> radii
        );

        //- Construct from dictionary
        cylinderListToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        cylinderListToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~cylinderListToCell();


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
