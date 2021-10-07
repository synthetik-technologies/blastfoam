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
    incompressibleMooneyRivlinElastic

Description
    Incompressible 3-parameter Mooney-Rivlin model, with strain-energy
    function defined as:

        psi = c10*(I1 - 3) + c01*(I2 - 3) + c11*(I1 - 3)*(I2 - 3)

    where I1 and I2 are the first and second invariants of the right
    Cauchy-Green deformation tensor; c10, c01, and c11 are material
    constants.

    The Cauchy stress tensor is given by:

        sigma = -p*I + 2*alpha*B + 2*beta*(B & B)

    where
    p        hydrostatic pressure (sigmaHyd == -p)
    I        indentity tensor
    B		left Cauchy-Green deformation tensor (== F & F.T)
    alpha, beta are coefficients that depend on material and the deformation:
            alpha =  c10 + c01*I1 + c11*(I1*(I1 - 3) + I2 -3)
            beta  = -c01 - c11*(I1 - 3)

    Incompressibility is enforced using the penalty method, where the bulk
    modulus (penalty) term should be set orders of magnitude greater than the
    shear modulus parameters.

    The bulk modulus acts as a penalty parameter, assuring that the
    incompressibility condition is satisfied. To ensure this, its value
    should be 'large enough' relative to the shear modulus parameters, but how
    large is subject to numerical experimentation. References below recommend it
    to be around 1e3 to 1e7 times the largest material parameter of the model
    (c01, c10, or c11). Based on numerical experiments, a value of 1e3 times the
    small-strain shear modulus of the Mooney-Rivlin model (mu = 2*(c01 + c10))
    should be enough for the current implementation.

    Optionally: a Poisson's equation can be solved to calculate sigmaHyd. This
    is required as incompressibility is approached.

    References
    - Costalat et al. (2011). Biomechanical wall properties of human
      intracranial aneurysms resected following surgical clipping
      (IRRAs Project). Journal of Biomechanics, 44(15), 2685–2691.
    - Nonlinear Solid Mechanics - A continuum approach for engineering
      Gerhard Hozapfel.


SourceFiles
    incompressibleMooneyRivlinElastic.C

Author
    Iago Lessa de Oliveira, FEIS/UNESP/UCD.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef incompressibleMooneyRivlinElastic_H
#define incompressibleMooneyRivlinElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class incompressibleMooneyRivlinElastic
:
    public mechanicalLaw
{
    // Private data

        // First material parameter
        const dimensionedScalar c10_;

        // First material parameter
        const dimensionedScalar c01_;

        // Third material parameter
        const dimensionedScalar c11_;

        // Shear modulus
        const dimensionedScalar mu_;

        // Bulk modulus
        const dimensionedScalar K_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        incompressibleMooneyRivlinElastic(const incompressibleMooneyRivlinElastic&);

        //- Disallow default bitwise assignment
        void operator=(const incompressibleMooneyRivlinElastic&);

public:

    //- Runtime type information
    TypeName("incompressibleMooneyRivlin");

    // Static data members


    // Constructors

        //- Construct from dictionary
        incompressibleMooneyRivlinElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~incompressibleMooneyRivlinElastic();


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
