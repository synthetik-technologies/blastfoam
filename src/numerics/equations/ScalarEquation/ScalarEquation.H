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
    Foam::scalarEquation

Description
    Base class for the a scalar equation

\*---------------------------------------------------------------------------*/

#ifndef ScalarEquation_H
#define ScalarEquation_H

#include "Equation.H"
#include "scalarEquation.H"
#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class scalarEquation Declaration
\*---------------------------------------------------------------------------*/

template<class InType>
class ScalarEquation
:
    public Equation<InType, scalar>,
    public scalarEquation
{
protected:
// Protected data

    //- Calculate the gradient using a finite difference
    void calculateGradient
    (
        const InType& x0,
        const label li,
        const scalar& fx0,
        scalarField& grad
    ) const;

    
public:

    typedef Equation<InType, scalar> equationType;
    typedef typename Equation<InType, scalar>::inType inType;
    typedef typename Equation<InType, scalar>::outType outType;
    
    // Constructors

        //- Construct given bounds
        ScalarEquation
        (
            const label nVar, 
            const InType& lowerLimit, 
            const InType& upperLimit
        );

        //- Construct given bounds
        ScalarEquation
        (
            const label nVar,
            const label nEqn, 
            const InType& lowerLimit, 
            const InType& upperLimit
        );


    //- Destructor
    virtual ~ScalarEquation();


    // Member Functions

        //- Return the number of input variables
        virtual label nVar() const
        {
            return equationType::nVar();
        }

        //- Return the number of equations
        virtual label nEqns() const 
        {
            return equationType::nEqns();
        }

        //- Return the number of derivatives
        virtual label nDerivatives() const = 0;

        //- Return the lower limit
        virtual tmp<scalarField> lowerLimits() const
        {
            return equationType::lowerLimits();
        }

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const
        {
            return equationType::upperLimits();
        }

        //- Return the finite difference step size for a direction
        virtual scalar dx(const label cmpti) const
        {
            return equationType::dx(cmpti);
        }
        
        //- Return the finite difference step size for a direction
        virtual scalar& dx(const label cmpti)
        {
            return equationType::dx(cmpti);
        }

        //- Return the finite difference step size
        virtual tmp<scalarField> dx() const
        {
            return equationType::dx();
        }

        //- Limit the value of x to the range of the equation
        virtual  void limit(InType& x) const
        {
            equationType::limit(x);
        }

        //- Calulate the gradient
        virtual void gradient
        (
            const InType& x,
            const label li,
            scalar& fx,
            scalarField& grad
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ScalarEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
