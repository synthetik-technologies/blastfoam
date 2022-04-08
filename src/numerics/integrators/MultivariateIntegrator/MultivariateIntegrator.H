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
    Foam::MultivariateIntegrator

Description
    Generic class for integrating a given univariate equation.

SourceFiles
    MultivariateIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef MultivariateIntegrator_H
#define MultivariateIntegrator_H

#include "multivariateIntegrator.H"
#include "univariateEquation.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class MultivariateIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MultivariateIntegrator
:
    public multivariateIntegrator
{
protected:

    // Protected data

        typedef univariateEquation<Type> equationType;

        //- Reference to Equation
        UautoPtr<const equationType> eqnPtr_;

        // Return the error give two values and spacing
        bool converged
        (
            const Type& t1,
            const Type& t2,
            const scalar dx,
            const scalar tol,
            const label i
        ) const
        {
            if ((dx - minDx_[i]) < 0)
            {
                return true;
            }
            else if (!adaptive())
            {
                return false;
            }
            return
                adaptiveError::error<Type>(t1, t2)
              < max(adaptiveError::mag<Type>(t1), absTolerance_[i])*tol;
        }

        virtual Type integrateFunc
        (
            const scalarField& x0,
            const scalarField& x1,
            const label li
        ) const = 0;

        //- Recursive add to the sum
        virtual void addMidsToInt
        (
            const label diri,
            label& fi,
            scalarField& x0,
            scalarField& x1,
            PtrList<Type>& fs,
            const label li
        ) const;

        //- Recursive add to the sum
        virtual void integrate_
        (
            const PtrList<Type>& Qs,
            const label diri,
            label& fi,
            scalarField& x0,
            scalarField& x1,
            const scalarField& tol,
            const label li,
            Type& fx
        ) const;


        //- Recursive integration
        virtual Type integrate_
        (
            const Type& Q,
            const label diri,
            const scalarField& x0,
            const scalarField& x1,
            const scalarField& tol,
            const label li
        ) const;

public:

    //- Runtime type information
    TypeName("MultivariateIntegrator");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            MultivariateIntegrator,
            dictionary,
            (const equationType& eqn, const dictionary& dict),
            (eqn, dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            MultivariateIntegrator,
            inputs,
            (
                const equationType& eqn,
                const multivariateIntegrator& inter
            ),
            (eqn, inter)
        );


    // Constructors

        //- Construct for given MultivariateIntegrator
        MultivariateIntegrator
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Copy constructor with a new equation
        MultivariateIntegrator
        (
            const equationType& eqn,
            const multivariateIntegrator& inter
        );

        //- Disallow default bitwise copy construction
        MultivariateIntegrator
        (
            const MultivariateIntegrator<Type>&
        ) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<MultivariateIntegrator<Type>> New
        (
            const equationType& eqn,
            const dictionary& dict
        );

        //- Select null constructed
        static autoPtr<MultivariateIntegrator<Type>> New
        (
            const equationType& eqn,
            const multivariateIntegrator& inter
        );


    //- Destructor
    virtual ~MultivariateIntegrator()
    {}


    // Member Functions

        //- Return if the integration is adaptive
        bool adaptive() const
        {
            return adaptiveError::adaptive<Type>(this->adaptive_);
        }

        //- Set the equation to integrate
        void setEqn(const equationType& eqn)
        {
            eqnPtr_.reset(&eqn);
        }

        //- Integrate the equation from x0 to x1
        virtual Type integrate
        (
            const scalarField& x0,
            const scalarField& x1,
            const label li = -1
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MultivariateIntegrator<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "MultivariateIntegrator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeMultivariateIntegrator(Type)                                       \
    defineNamedTemplateTypeNameAndDebug(MultivariateIntegrator<Type>, 0);      \
                                                                               \
    defineTemplateRunTimeSelectionTable(MultivariateIntegrator<Type>, dictionary);\
    defineTemplateRunTimeSelectionTable(MultivariateIntegrator<Type>, inputs);


#define makeMultivariateIntegratorType(MultivariateIntegratorType, Type)       \
    defineNamedTemplateTypeNameAndDebug(MultivariateIntegratorType<Type>, 0);  \
                                                                               \
    MultivariateIntegrator<Type>::                                             \
        adddictionaryConstructorToTable<MultivariateIntegratorType<Type>>      \
        add##MultivariateIntegratorType##Type##dictionary##ConstructorToTable_;\
                                                                               \
    MultivariateIntegrator<Type>::                                             \
        addinputsConstructorToTable<MultivariateIntegratorType<Type>>          \
        add##MultivariateIntegratorType##Type##inputs##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
