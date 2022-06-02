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
    Foam::CoefficientUnivariateEquation

Description
    Base class for holding modifiable coefficents

\*---------------------------------------------------------------------------*/

#ifndef Coefficients_H
#define Coefficients_H

#include "UnivariateEquation.H"
#include "Equation.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class Coefficients Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Coefficients
{
protected:
    // Protected member data

        //- List of coefficients
        Field<scalar> coeffs_;


public:

    typedef UList<scalar> XType;

    // Constructors

        //- Construct given number of coeffients
        Coefficients(const label nCoeffs)
        :
            coeffs_(nCoeffs)
        {}

        //- Construct given coeffients
        Coefficients(const List<scalar>& coeffs)
        :
            coeffs_(coeffs)
        {}


    //- Destructor
    virtual ~Coefficients()
    {}

    // Member functions

        //- Set the number of coefficients
        inline void setN(const label n)
        {
            coeffs_.setSize(n);
        }

        //- Return the number of coefficients
        inline label n() const
        {
            return coeffs_.size();
        }

        //- Return the coefficients
        const Field<scalar>& coeffs() const
        {
            return coeffs_;
        }

        //- Access the coefficients
        Field<scalar>& coeffsRef()
        {
            return coeffs_;
        }

        //- Return the jacobian w.r.t. the coefficients
        virtual void coeffJ
        (
            const XType& x,
            const label i,
            RectangularMatrix<Type>& J,
            const label li
        ) const
        {
            NotImplemented;
        }
};


/*---------------------------------------------------------------------------*\
                        Class CoefficientEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CoefficientEquation
:
    public Equation<Type>,
    public Coefficients<Type>
{
public:

    // Constructors

        //- Construct given number of equations
        CoefficientEquation
        (
            const label nCoeffs,
            const string& eqnString = string::null
        )
        :
            Equation<Type>(eqnString),
            Coefficients<Type>(nCoeffs)
        {}

        //- Construct given bounds
        CoefficientEquation
        (
            const scalar& lowerLimit,
            const scalar& upperLimit,
            const List<scalar>& coeffs,
            const string& eqnString = string::null
        )
        :
            Equation<Type>(lowerLimit, upperLimit, eqnString),
            Coefficients<Type>(coeffs)
        {}

        //- Construct from dictionary
        CoefficientEquation(const dictionary& dict, const label nCoeffs)
        :
            UnivariateEquation<Type>(dict),
            Coefficients<Type>(dict)
        {}


    //- Destructor
    virtual ~CoefficientEquation()
    {}
};


template<class Type>
class LinearEquation
:
    public CoefficientEquation<Type>
{
public:

    //- Constructor

        //- Construct given number of equations
        LinearEquation
        (
            const string& eqnString = string::null
        )
        :
            CoefficientEquation<Type>(2, eqnString)
        {}

        //- Construct given bounds
        LinearEquation
        (
            const scalar& lowerLimit,
            const scalar& upperLimit,
            const string& eqnString = string::null
        )
        :
            CoefficientEquation<Type>
            (
                lowerLimit,
                upperLimit,
                2,
                eqnString
            )
        {}

        //- Construct from dictionary
        LinearEquation(const dictionary& dict)
        :
            CoefficientEquation<Type>(dict)
        {
            this->setN(2);
        }

    // Member functions

        label nDerivatives() const
        {
            return 1;
        }

        scalar fx(const scalar x, const label li) const
        {
            return this->coeffs_[0] + this->coeffs_[1]*x;
        }

        //- Calculate the gradient
        virtual scalar dfdx(const scalar x, const label li) const
        {
            return this->coeffs_[1];
        }
};



/*---------------------------------------------------------------------------*\
                Class UnivariateCoefficientEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UnivariateCoefficientEquation
:
    public UnivariateEquation<Type>,
    public Coefficients<Type>
{
public:

    // Constructors

        //- Construct given number of equations
        UnivariateCoefficientEquation
        (
            const label nVar,
            const label nCoeffs,
            const string& eqnString = string::null
        )
        :
            UnivariateEquation<Type>(nVar, eqnString),
            Coefficients<Type>(nCoeffs)
        {}

        //- Construct given bounds
        UnivariateCoefficientEquation
        (
            const scalarList& lowerLimit,
            const scalarList& upperLimit,
            const List<scalar>& coeffs,
            const string& eqnString = string::null
        )
        :
            UnivariateEquation<Type>(lowerLimit, upperLimit, eqnString),
            Coefficients<Type>(coeffs)
        {}

        //- Construct from dictionary
        UnivariateCoefficientEquation(const dictionary& dict, const label nCoeffs)
        :
            UnivariateEquation<Type>(dict),
            Coefficients<Type>(nCoeffs)
        {}


    //- Destructor
    virtual ~UnivariateCoefficientEquation()
    {}
};


template<class Type>
class LinearUnivariateEquation
:
    public UnivariateCoefficientEquation<Type>
{
public:

    //- Constructor

        //- Construct given number of equations
        LinearUnivariateEquation
        (
            const label nVar,
            const string& eqnString = string::null
        )
        :
            UnivariateCoefficientEquation<Type>(nVar, nVar + 1, eqnString)
        {}

        //- Construct given bounds
        LinearUnivariateEquation
        (
            const scalarList& lowerLimit,
            const scalarList& upperLimit,
            const string& eqnString = string::null
        )
        :
            UnivariateCoefficientEquation<Type>
            (
                lowerLimit,
                upperLimit,
                lowerLimit.size() + 1,
                eqnString
            )
        {}

        //- Construct from dictionary
        LinearUnivariateEquation(const dictionary& dict)
        :
            UnivariateCoefficientEquation<Type>(dict)
        {
            this->setN(this->nVar_ + 1);
        }

    // Member functions

        label nDerivatives() const
        {
            return 1;
        }

        scalar fX
        (
            const typename UnivariateCoefficientEquation<Type>::VarType& x,
            const label li
        ) const
        {
            scalar res = this->coeffs_[0];
            forAll(x, i)
            {
                res += this->coeffs_[i+1]*x[i];
            }
            return res;
        }

        //- Calculate the gradient
        virtual void dfdX
        (
            const typename UnivariateCoefficientEquation<Type>::VarType& x,
            const label li,
            List<Type>& dfdx
        ) const
        {
            dfdx.setSize(x.size());
            forAll(x, i)
            {
                dfdx[i] = this->coeffs_[i+1];
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
