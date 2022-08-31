/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    cohesiveFvPatch

Description
    Cohesive fv patch.

SourceFiles
    cohesiveFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cohesiveFvPatch_H
#define cohesiveFvPatch_H

#include "fvPatch.H"
#include "surfaceFields.H"
#include "cyclicLduInterface.H"
#include "cohesivePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cohesiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cohesiveFvPatch
:
    public fvPatch
{
    // Private data

        const cohesivePolyPatch& cohesivePolyPatch_;


protected:

    // Protected Member Functions

        //- Make patch face non-orthogonality correction vectors
        virtual void makeCorrVecs(fvsPatchVectorField&) const;


public:

    //- Runtime type information
    TypeName(cohesivePolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cohesiveFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            fvPatch(patch, bm),
            cohesivePolyPatch_
            (
                refCast<const cohesivePolyPatch>(patch)
            )
        {}


    // Destructor

        virtual ~cohesiveFvPatch()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
