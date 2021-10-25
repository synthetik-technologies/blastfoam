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
    Abstract base-class for a integrator

SourceFiles
    Integrator.C

\*---------------------------------------------------------------------------*/

#ifndef Integrator_H
#define Integrator_H

#include "univariateEquation.H"
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
{
protected:

    // Protected data

        typedef univariateEquation<Type> equationType;

        //- Reference to Equation
        UautoPtr<const equationType> eqnPtr_;

        //- Number of steps to integrate with
        label nSteps_;

        //- Number of intervals
        label nIntervals_;

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


    // Constructors

        //- Construct for given Integrator
        Integrator(const equationType& eqn, const dictionary& dict);

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
            const word& intType
        );


    //- Destructor
    virtual ~Integrator()
    {}


    // Member Functions

        //- Return the number of steps
        label nSteps() const
        {
            return nSteps_;
        }

        //- Return the number of intervals
        label nIntervals() const
        {
            return nIntervals_;
        }

        //- Access the number of steps
        virtual void setNIntervals(const label nIntervals)
        {
            nSteps_ = nIntervals;
            nIntervals_ = nIntervals;
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
    defineTemplateRunTimeSelectionTable(Integrator<Type>, dictionary);


#define makeIntegratorType(IntegratorType, Type)                               \
    defineNamedTemplateTypeNameAndDebug(IntegratorType<Type>, 0);              \
                                                                               \
    Integrator<Type>::adddictionaryConstructorToTable<IntegratorType<Type>>    \
    addIntegrator##IntegratorType##Type##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
