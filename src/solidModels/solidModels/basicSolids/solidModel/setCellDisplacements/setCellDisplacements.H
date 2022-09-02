/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    setCellDisplacements

Description
    List of displacement to be enforced as given cell locations, via the
    fvMatrix "setValues" function.

    In this way, displacements can be enforced a locations internal to the mesh
    rather than just at bondary conditions.

    The displacement for each cell can be constant in time or time-varying.

Author
    Philip Cardiff, UCD.

SourceFiles
    setCellDisplacements.C

\*---------------------------------------------------------------------------*/

#ifndef setCellDisplacements_H
#define setCellDisplacements_H

#include "dictionary.H"
#include "fvMesh.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class setCellDisplacements Declaration
\*---------------------------------------------------------------------------*/

class setCellDisplacements
{
    // Private data

        //- Const-reference to the mesh
        const fvMesh& mesh_;

        //- Copy of the dictionary
        const dictionary dict_;

        //- List of cell indices to be set
        labelList cellIDs_;

        //- Constant displacements for the cells
        vectorField constantDisps_;

        //- Time-varying displacements for the cells
        PtrList<Function1s::Table<vector> > timeVaryingDisps_;

        //- Current cell displacements
        mutable vectorField currentCellDisps_;

        //- Current time-step index
        mutable label curTimeIndex_;


    // Private Member Functions

        //- Read settings from the dictionary and initialise fields
        void readDict();

        //- Disallow default bitwise copy construct
        setCellDisplacements(const setCellDisplacements&);

        //- Disallow default bitwise assignment
        void operator=(const setCellDisplacements&);

public:

    //- Runtime type information
    TypeName("setCellDisplacements");


    // Constructors

        //- Construct from dictionary
        setCellDisplacements
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~setCellDisplacements()
        {}


    // Member Functions

        // Evaluation

            //- Return a const-reference to cellIDs
            const labelList& cellIDs() const
            {
                return cellIDs_;
            }

            //- Return a const-reference to the cell displacements for the
            //  current time
            const vectorField& cellDisps() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
