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
    viscousHookeanElastic

Description
    Visco-elastic Hookean mechanical law.

    The model formulation and stress calculation follow the procedure given in:

        Simo and Hughes, 2000, Computational Inelasticity, Springer

    The material is represented as a generalised Maxwell model (aka Prony
    series) i.e. a spring in parallel with multiple Maxwell models, where a
    Maxwell model is a spring in series with a dashpot.

    The relaxation times can be optionally scaled by the temperature field using
    the Williams-Landel-Ferry (WLF) approximation and the parameters C1 and C2,
    where the scaling factor is: a(T) = -C1*(T - T0)/(C2 + (T - T0));

    The material is specifed as, for example:

    @verbatim
    plastic
    {
        type            viscousHookeanElastic;
        EInfinity       EInfinity [1 -1 -2 0 0 0 0] 200e+9;
        E               ( 100e+09 50e+09 10e+09 );
        relaxationTimes ( 100 10 1 );
        nu              nu [0 0 0 0 0 0 0] 0.3;
        //WilliamsLandelFerry yes;
        //C1            C1 [0 0 0 0 0 0 0] 17.44;
        //C2            C2 [0 0 0 1 0 0 0] 51.6;
    }
    @endverbatim

SourceFiles
    viscousHookeanElastic.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef viscousHookeanElastic_H
#define viscousHookeanElastic_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class viscousHookeanElastic Declaration
\*---------------------------------------------------------------------------*/

class viscousHookeanElastic
:
    public mechanicalLaw
{
    // Private data

        //- Relaxed (long term) modulus of elasticity
        dimensionedScalar EInf_;

        //- List of modulii of elasticity for each Maxwell model
        List<scalar> E_;

        //- List of relaxation times for each Maxwell model
        //  Note: tau_i = eta_i/E_i, where eta_i is the viscosity of the dashpot
        //  for the i^th Maxwell model
        List<scalar> tau_;

        //- Relaxed (long term) nondimensional relative modulus
        dimensionedScalar gammaInf_;

        //- List of nondimensional relative modulii
        List<scalar> gamma_;

        //- Relaxed (long term) Poisson's ratio
        dimensionedScalar nu_;

        //- Relaxed (long term) first Lame parameter
        dimensionedScalar lambda_;

        //- Relaxed (long term) second Lame parameter (shear modulus)
        dimensionedScalar mu_;

        //- Relaxed (long term) and unrelaxed (short term) bulk modulus
        //  Only the deviatoric component of stress is viscous
        dimensionedScalar k_;

        //- Internal stress variables for each Maxwell model
        PtrList<volSymmTensorField> h_;

        //- Internal stress variables for each Maxwell model surface fields
        PtrList<surfaceSymmTensorField> hf_;

        //- Deviatoric component of the "initial stress"
        volSymmTensorField s_;

        //- Deviatoric component of the "initial stress" surface field
        surfaceSymmTensorField sf_;

        //- Enable WLF method to scale the relaxation times based on the
        //  temperature field
        const Switch WilliamsLandelFerryShift_;

        //- First WLF parameter
        const dimensionedScalar C1_;

        //- Second WLF parameter
        const dimensionedScalar C2_;

        //- Reference temperature for WLF
        const dimensionedScalar Tref_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        viscousHookeanElastic(const viscousHookeanElastic&);

        //- Disallow default bitwise assignment
        void operator=(const viscousHookeanElastic&);

public:

    //- Runtime type information
    TypeName("viscousHookeanElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        viscousHookeanElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~viscousHookeanElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return the implicit stiffness for a patch
        //  This is the diffusivity for the Laplacian term
        virtual tmp<scalarField> impK(const label) const;

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return elastic modulus
        virtual tmp<volScalarField> elasticModulus() const;

        //- Return shear modulus
        virtual tmp<volScalarField> shearModulus() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Return the desired new time-step
        //virtual scalar newDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
