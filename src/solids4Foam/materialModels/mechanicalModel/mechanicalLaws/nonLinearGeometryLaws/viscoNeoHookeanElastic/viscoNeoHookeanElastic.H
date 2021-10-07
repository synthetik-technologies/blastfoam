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
    viscoNeoHookeanElastic

Description
    Neo-visco-elastic Hookean mechanical law:

        tau_ = J*p*I + g(deltaT)*tau_+hBar

    The model formulation and stress calculation follow the procedure given in:

        Simo and Hughes, 2000, Computational Inelasticity, Springer

    The material is represented as a generalised Maxwell model (aka Prony
    series) i.e. a spring in parallel with multiple Maxwell models, where a
    Maxwell model is a spring in series with a dashpot.

SourceFiles
    viscoNeoHookeanElastic.C

Author
    Philip Cardiff, UCD. All rights reserved
    Xiaoxue Shen, CSU.

\*---------------------------------------------------------------------------*/

#ifndef viscoNeoHookeanElastic_H
#define viscoNeoHookeanElastic_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class viscoNeoHookeanElastic Declaration
\*---------------------------------------------------------------------------*/

class viscoNeoHookeanElastic
:
    public mechanicalLaw
{
    // Private data

        //- Relaxed (long term) modulus of elasticity
        dimensionedScalar EInf_;

        //- List of modulii of elasticity for each Maxwell model
        List<scalar> E_;

        //- List of relaxation times for each Maxwell model
        List<scalar> tau_;

        //- Relaxed (long term) nondimensional relative modulus
        scalar gammaInf_;

        //- List of nondimensional relative modulii
        List<scalar> gamma_;

        //- Relaxed (long term) Poisson's ratio
        dimensionedScalar nu_;

        //- Relaxed (long term) first Lame parameter
        dimensionedScalar lambda_;

        //- Relaxed (long term) second Lame parameter (shear modulus)
        dimensionedScalar mu_;

        //- Shear modulii in Ogden model
        dimensionedScalar mu0_;
        dimensionedScalar mu1_;
        dimensionedScalar mu2_;

        //- Dimensionless material constants
        List<scalar> alpha_;

        //- An empirical coefficient
        dimensionedScalar beta_;

        //- Relaxed (long term) and unrelaxed (short term) bulk modulus
        //  Only the deviatoric component of stress is viscous
        dimensionedScalar K_;


        //- Internal stress variables for each Maxwell model
        PtrList<volSymmTensorField> h_;

        //- Internal stress variables for each Maxwell model surface fields
        PtrList<surfaceSymmTensorField> hf_;

        //- Internal stress variables (intermediate variables,
        //  H_=exp(-deltaT/tau)*h_-exp(-deltaT/(2*tau))*s_) for each Maxwell
        //  model
        PtrList<volSymmTensorField> H_;

        //- Internal stress variables, intermediate variables,
        //  for each Maxwell model surface fields
        PtrList<surfaceSymmTensorField> Hf_;

        //- Intermediate variables, won't be stored
        PtrList<volSymmTensorField> transformH_;
        PtrList<surfaceSymmTensorField> transformHf_;

        //- Deviatoric component of the "initial stress"
        volSymmTensorField s_;

        //- Deviatoric component of the "initial stress" surface field
        surfaceSymmTensorField sf_;

        //- Intermediate variables, won't be stored
        volSymmTensorField transformNeeded_;
        volSymmTensorField transformFbar_;

        //- Intermediate variables, won't be stored
        surfaceSymmTensorField transformNeededf_;
        surfaceSymmTensorField transformFbarf_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        viscoNeoHookeanElastic(const viscoNeoHookeanElastic&);

        //- Disallow default bitwise assignment
        void operator=(const viscoNeoHookeanElastic&);

public:

    //- Runtime type information
    TypeName("viscoNeoHookeanElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        viscoNeoHookeanElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~viscoNeoHookeanElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return the bulk modulus
        virtual tmp<volScalarField> K() const;

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
