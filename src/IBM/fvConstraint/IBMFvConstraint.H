/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::fv::IBMForceConstraint

Description
    Adds immersed boundary forcing

Usage
    Example usage:
    \verbatim
    ibm
    {
        type        IBM;
    }
    \endverbatim

SourceFiles
    IBMFvConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef IBMForceConstraint_H
#define IBMForceConstraint_H

#include "fvModel.H"
#include "fvConstraint.H"
#include "immersedBoundaryObjectListSolver.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class IBMForceConstraint Declaration
\*---------------------------------------------------------------------------*/

class IBMForceConstraint
:
    public fvConstraint
{
    // Private Data

        //- Immersed objects
        immersedBoundaryObjectListSolver& ibm_;

        //- Fields with IBM forcing added
        hashedWordList forcedFields_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the phase density (alpha*rho)
        tmp<volScalarField> alphaRho(const word& phase) const;

        //- Return the phase old density (alpha*rho)
        tmp<volScalarField> alphaRhoOld(const word& phase) const;

        //- Add source terms to a phase equation
        template<class Type>
        bool constrainType
        (
            fvMatrix<Type>& eqn,
            const word& fieldName
        ) const;


public:

    //- Runtime type information
    TypeName("IBM");


    // Constructors

        //- Construct from explicit source name and mesh
        IBMForceConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        IBMForceConstraint(const IBMForceConstraint&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList constrainedFields() const;


        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_CONSTRAINT_CONSTRAIN);

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const IBMForceConstraint&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IBMFvConstraintTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
