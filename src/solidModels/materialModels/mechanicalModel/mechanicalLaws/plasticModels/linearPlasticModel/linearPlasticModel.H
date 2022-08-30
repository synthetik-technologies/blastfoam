/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
31-03-2022 Synthetik Applied Technologies: Added base plasticModel
-------------------------------------------------------------------------------
License
    This file is a derivative work of foam-extend.

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
    linearPlasticModel

Description
    Base class for a plastic model with linear elastic strain

SourceFiles
    linearPlasticModel.C

\*---------------------------------------------------------------------------*/

#ifndef linearPlasticModel_H
#define linearPlasticModel_H

#include "mechanicalLaw.H"
#include "plasticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class linearPlasticModel
:
    public plasticModel
{
protected:


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        linearPlasticModel(const linearPlasticModel&);

        //- Disallow default bitwise assignment
        void operator=(const linearPlasticModel&);

        //- Calculate hydrostatic component of the stress tensor
        tmp<volScalarField> hydrostaticStress
        (
            const volScalarField& trEpsilon
        );

        //- Calculate hydrostatic component of the stress tensor
        tmp<surfaceScalarField> hydrostaticStress
        (
            const surfaceScalarField& trEpsilon
        );


public:

    //- Runtime type information
    TypeName("linearPlasticModel");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearPlasticModel
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~linearPlasticModel();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return the implicit stiffness for a patch
        //  This is the diffusivity for the Laplacian term
        virtual tmp<scalarField> impK(const label) const;

        //- Update the stress
        virtual void correct(volSymmTensorField& sigma)
        {
            correct(sigma, true);
        }

        virtual void correct
        (
            volSymmTensorField& sigma,
            const bool needUpdate
        );

        //- Update the stress surface field
        virtual void correct(surfaceSymmTensorField& sigma)
        {
            correct(sigma, true);
        }

        virtual void correct
        (
            surfaceSymmTensorField& sigma,
            const bool needUpdate
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
