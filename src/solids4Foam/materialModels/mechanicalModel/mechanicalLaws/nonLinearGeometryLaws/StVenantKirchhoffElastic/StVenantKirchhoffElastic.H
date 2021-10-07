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
    StVenantKirchhoffElastic

Description
    St. Venant Kirchhoff Hyperelastic law, where the 2nd Piola Kirchhoff stress
    is given as:

        S = 2*mu*G + lambda*tr(G)*I

    where
    G        Green strain
    lambda   first Lame parameter
    mu       second Lame parameter (shear modulus)

    Either E (Young's modulus) and nu (Poisson's ratio) can be provided or mu
    (shear modulus) and K (bulk modulus)

SourceFiles
    StVenantKirchhoffElastic.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef StVenantKirchhoffElastic_H
#define StVenantKirchhoffElastic_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class StVenantKirchhoffElastic
:
    public mechanicalLaw
{
    // Private data

        // First Lame parameter
        dimensionedScalar lambda_;

        // Second Lame parameter (shear modulus)
        dimensionedScalar mu_;

        // Bulk modulus
        dimensionedScalar K_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        StVenantKirchhoffElastic(const StVenantKirchhoffElastic&);

        //- Disallow default bitwise assignment
        void operator=(const StVenantKirchhoffElastic&);

public:

    //- Runtime type information
    TypeName("StVenantKirchhoffElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        StVenantKirchhoffElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~StVenantKirchhoffElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

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
