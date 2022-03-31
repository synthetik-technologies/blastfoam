/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    poroLinearElastic

Description
    Poro-elastic mechanical law where the total stress is the sum of the solid
    skeleton (effective) stress, governed by Hooke's law, and the poro-pressure.

    More reading at:

    T.Tang, O.Hededal & P.Cardiff (2014). On finite volume method implementation
    of poro-elasto-plasticity soil model. International Journal for
    Numerical and Analytical Methods in Geomechanics, 10.1002/nag.2361.

    T.Tang & O.Hededal (2014). Simulation of pore pressure accumulation under
    cyclic loading using finite volume method. Proceedings of the 8th European
    Conference on Numerical Methods in Geotechnical Engineering (NUMGE14),
    Volume 2, Pages 1301-1306.

SourceFiles
    poroLinearElastic.C

Author
    Tian Tang, DTU. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef poroLinearElastic_H
#define poroLinearElastic_H

#include "linearElastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poroLinearElastic Declaration
\*---------------------------------------------------------------------------*/

class poroLinearElastic
:
    public linearElastic
{
    // Private Data

        //- Initial (residual) poro-pressure field
        //  Defaults to zero if not specified
        const volScalarField p0_;

        //- Initial (residual) stress surface field
        mutable surfaceScalarField* p0fPtr_;

    // Private Member Functions

        //- Make the initial pore-pressure surface field
        void makeP0f() const;

        //- Return a reference to the initial pore-pressure surface field
        const surfaceScalarField& p0f() const;

        //- Disallow default bitwise copy construct
        poroLinearElastic(const poroLinearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const poroLinearElastic&);

public:

    //- Runtime type information
    TypeName("poroLinearElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        poroLinearElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~poroLinearElastic();


    // Member Functions

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
