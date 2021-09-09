/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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
    Foam::massToCell

Description
    Basic class used for mass based cell sets

SourceFiles
    massToCell.C

\*---------------------------------------------------------------------------*/

#ifndef massToCell_H
#define massToCell_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class massToCell Declaration
\*---------------------------------------------------------------------------*/

class massToCell
{

    // Private Data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Should the dictionary be read
        bool mustRead_;

        //- Density field
        scalar rho_;

        //- Mass to set
        scalar mass_;

        //- Centre of mass
        vector centre_;

        //- Volume to set
        scalar volume_;

protected:

    // has the dictionary been read
    mutable bool read_;

public:

    // Constructors

        //- Construct from dictionary
        massToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~massToCell();

    //- Do values need to be read
    bool mustRead() const
    {
        return mustRead_;
    }

    //- Return the volume needed
    scalar volume() const
    {
        return volume_;
    }

    //- Return the center of mass
    vector centre() const
    {
        return centre_;
    }

    //- Return the center of mass
    vector& centre()
    {
        return centre_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
