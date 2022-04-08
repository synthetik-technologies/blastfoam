/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::Integrator

Description
    Generic class for integrating a given univariate equation.

SourceFiles
    Integrator.C

\*---------------------------------------------------------------------------*/

#ifndef Integrator_H
#define Integrator_H

#include "integrator.H"
#include "Equation.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Integrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Integrator
:
    public integrator
{
protected:

    // Protected data

        typedef Equation<Type> equationType;

        //- Reference to Equation
        UautoPtr<const equationType> eqnPtr_;

        // Return the error give two values and spacing
        bool converged
        (
            const Type& t1,
            const Type& t2,
            const scalar dx,
            const scalar tol
        ) const
        {
            if (mag(dx) < minDx_)
            {
                return true;
            }
            else if (!adaptive())
            {
                return false;
            }
            return
                adaptiveError::error<Type>(t1, t2)
              < max(adaptiveError::mag<Type>(t1), absTolerance_)*tol;
        }


public:

    //- Runtime type information
    TypeName("Integrator");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            Integrator,
            dictionary,
            (const equationType& eqn, const dictionary& dict),
            (eqn, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            Integrator,
            inputs,
            (
                const equationType& eqn,
                const integrator& inter
            ),
            (eqn, inter)
        );


    // Constructors

        //- Construct for given Integrator
        Integrator(const equationType& eqn, const dictionary& dict);

        //- Copy constructor with a new equation
        Integrator
        (
            const equationType& eqn,
            const integrator& inter
        );

        //- Disallow default bitwise copy construction
        Integrator(const Integrator<Type>&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<Integrator<Type>> New
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Select null constructed
        static autoPtr<Integrator<Type>> New
        (
            const equationType& eqn,
            const integrator& inter
        );


    //- Destructor
    virtual ~Integrator()
    {}


    // Member Functions

        //- Return if the integration is adaptive
        inline bool adaptive() const
        {
            return adaptiveError::adaptive<Type>(adaptive_);
        }

        //- Set the equation to integrate
        void setEqn(const equationType& eqn)
        {
            eqnPtr_.reset(&eqn);
        }

        //- Integrate the equation from x0 to x1
        virtual Type integrate
        (
            const scalar x0,
            const scalar x1,
            const label li = -1
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Integrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "Integrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIntegrator(Type)                                                   \
    defineNamedTemplateTypeNameAndDebug(Integrator<Type>, 0);                  \
                                                                               \
    defineTemplateRunTimeSelectionTable(Integrator<Type>, dictionary);         \
    defineTemplateRunTimeSelectionTable(Integrator<Type>, inputs);


#define makeIntegratorType(IntegratorType, Type)                               \
    defineNamedTemplateTypeNameAndDebug(IntegratorType<Type>, 0);              \
                                                                               \
    Integrator<Type>::adddictionaryConstructorToTable<IntegratorType<Type>>    \
    add##IntegratorType##Type##dictinoary##ConstructorToTable_;                \
                                                                               \
    Integrator<Type>::addinputsConstructorToTable<IntegratorType<Type>>        \
    add##IntegratorType##Type##inputs##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
