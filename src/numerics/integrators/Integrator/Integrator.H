/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "integratorBase.H"
#include "Equation.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class integrator Declaration
\*---------------------------------------------------------------------------*/

class integrator
:
    public integratorBase
{
protected:
    // Protected data

        //- Tolerance
        scalar tolerance_;

        //- Absolute tolerance
        scalar absTolerance_;

        //- Maximum number of domain splits
        label maxSplits_;

        //- Number of intervals
        label nIntervals_;

        mutable label intervals_;
        mutable scalar minDx_;

        //- Reset counters and minimum spacing
        void reset(const scalar dx) const;


public:

    //- Runtime type information
    TypeName("integrator");

    // Constructors

        //- Construct for from dictinoary
        integrator(const dictionary& dict);


    //- Destructor
    virtual ~integrator()
    {}


    // Member Functions

        //- Return if the integration is adaptive
        virtual bool adaptive() const = 0;

        //- Return the number of intervals
        label nIntervals() const
        {
            return adaptive() ? intervals_ : nIntervals_;
        }

        //- Access the number of intervals to use
        virtual void setNIntervals(const label nIntervals)
        {
            adaptive_ = false;
            nIntervals_ = nIntervals;
        }
};

/*---------------------------------------------------------------------------*\
                          Class Integrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Adapt>
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
                adaptiveError::error<Adapt>
                (
                    static_cast<const Adapt&>(t1),
                    static_cast<const Adapt&>(t2)
                )
              < max
                (
                    adaptiveError::mag<Adapt>(static_cast<const Adapt&>(t1)),
                    absTolerance_
                )*tol;
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
            (
                const equationType& eqn,
                const dictionary& dict
            ),
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
        Integrator(const Integrator<Type, Adapt>&) = delete;


    // Selectors

        //- Construct from dictionary
        static autoPtr<Integrator<Type, Adapt>> New
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Construct from dictionary with a given type
        static autoPtr<Integrator<Type, Adapt>> New
        (
            const word& integratorTypeName,
            const equationType& eqn,
            const dictionary& dict
        );

        //- Select null constructed
        static autoPtr<Integrator<Type, Adapt>> New
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
            return adaptiveError::adaptive<Adapt>(adaptive_);
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
        void operator=(const Integrator<Type, Adapt>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IntegratorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIntegrator(Type, Adapt)                                            \
    typedef Integrator<Type, Adapt> Integrator##Type##Adapt;                   \
    defineNamedTemplateTypeNameAndDebug(Integrator##Type##Adapt, 0);           \
                                                                               \
    defineTemplateRunTimeSelectionTable(Integrator##Type##Adapt, dictionary);  \
                                                                               \
    defineTemplateRunTimeSelectionTable(Integrator##Type##Adapt, inputs);


#define makeIntegratorType(IntegratorType, Type, Adapt)                        \
    typedef IntegratorType<Type, Adapt> IntegratorType##Type##Adapt;           \
    defineNamedTemplateTypeNameAndDebug(IntegratorType##Type##Adapt, 0);       \
                                                                               \
    Integrator##Type##Adapt::adddictionaryConstructorToTable                   \
    <                                                                          \
        IntegratorType##Type##Adapt                                            \
    > add##IntegratorType##Type##Adapt##dictinoary##ConstructorToTable_;       \
                                                                               \
    Integrator##Type##Adapt::addinputsConstructorToTable                       \
    <                                                                          \
        IntegratorType##Type##Adapt                                            \
    > add##IntegratorType##Type##Adapt##inputs##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
