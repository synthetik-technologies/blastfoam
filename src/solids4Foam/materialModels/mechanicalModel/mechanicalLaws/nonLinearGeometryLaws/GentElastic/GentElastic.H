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
    GentElastic

Description
    Gent-Flory elasticity for the prediction of microgels behavior, where
    J        Jacobian of the deformation gradient
    I        indentity tensor
    F        deformation gradient

    References

    Flory & Huggins, Development of Expert Networks: A Hybrid System of Expert Systems and Neural Networks,
    Neural Networks in Bioprocessing and Chemical Engineering,1995.

    Gent, A new constitutive relation for rubber, Rubber Chemistry Technology, 1996.

SourceFiles
    GentElastic.C

Author
    Thomas Nevolianis, RWTH Aachen
    Some re-factoring by Philip Cardiff, UCD

\*---------------------------------------------------------------------------*/

#ifndef GentElastic_H
#define GentElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class GentElastic
:
    public mechanicalLaw
{
    // Private data

        // Young's modulus
        const dimensionedScalar E_;

        // Poisson's ratio
        const dimensionedScalar nu_;

        // Shear modulus
        const dimensionedScalar mu_;

        // Temperature
        const dimensionedScalar temperature_;

        // ilim = 3*( pow(Jlim,2/3) - 1 ) where Jlim is the maximum volume fraction
        const dimensionedScalar ilim_;

        // Number of crosslinks
        const dimensionedScalar N_;

        // Unit volume of solvent
        const dimensionedScalar Vs_;

        // Bulk modulus
        const dimensionedScalar K_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GentElastic(const GentElastic&);

        //- Disallow default bitwise assignment
        void operator=(const GentElastic&);

public:

    //- Runtime type information
    TypeName("GentElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        GentElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~GentElastic();


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
