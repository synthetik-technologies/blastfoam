/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
12-04-2022 Synthetik Applied Technologies : Added equation functionality
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
    Foam::CodedEquation

Description
    User defined equation with a single input and single output

\*---------------------------------------------------------------------------*/

#ifndef CodedCoeffientEquation_H
#define CodedCoeffientEquation_H

#include "Coefficients.H"
#include "dictionary.H"

#include "codedBase.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "OSspecific.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class CodedEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CodedCoefficientEquation
:
    public UnivariateCoefficientEquation<Type>,
    public codedBase
{
    // Private Data

        //- Number of derivatives
        label nDerivatives_;

        //- The dynamically generated Equation pointer
        mutable autoPtr<UnivariateCoefficientEquation<Type>> redirectEquationPtr_;


    // Private Member Functions

        //- Get the keywords associated with source code
        virtual wordList codeKeys() const;

         //- Get the name of the dictionary variables in the source code
        virtual wordList codeDictVars() const;

        //- Adapt the context for the current object
        virtual void prepare
        (
            dynamicCode& dynCode,
            const dynamicCodeContext& context
        ) const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const;

        //- Compile, link and return the now coded Equation
        autoPtr<UnivariateCoefficientEquation<Type>> compileNew();

        //- Expand strings in the code entry
        const dictionary& expandCodeDict(const dictionary& cDict) const;


public:

    TypeName("codedCoefficient");

    // Constructors

        //- Construct from dictionary
        CodedCoefficientEquation(const dictionary& dict);

        //- Construct from dictionary and object registry
        CodedCoefficientEquation(const objectRegistry& obr, const dictionary& dict);


    //- Destructor
    virtual ~CodedCoefficientEquation();


    // Member Functions

        //- Set the number of coefficients
        virtual inline void setN(const label n)
        {
            this->coeffs_.setSize(n);
            redirectEquationPtr_->setN(n);
        }

        //- Return the number of coefficients
        virtual inline label n() const
        {
            return redirectEquationPtr_->n();
        }

        //- Return the coefficients
        virtual const Field<scalar>& coeffs() const
        {
            return redirectEquationPtr_->coeffs();
        }

        //- Access the coefficients
        virtual Field<scalar>& coeffsRef()
        {
            return redirectEquationPtr_->coeffsRef();
        }

        //- Return the number of implemented derivatives
        virtual inline label nDerivatives() const
        {
            return redirectEquationPtr_->nDerivatives();
        }

        //- Return the function value
        virtual inline Type fX
        (
            const typename UnivariateEquation<Type>::VarType& x,
            const label li
        ) const
        {
            return redirectEquationPtr_->fX(x, li);
        }

        //- Return the function value
        virtual void dfdX
        (
            const typename UnivariateEquation<Type>::VarType& x,
            const label li,
            List<Type>& dfdx
        ) const
        {
            redirectEquationPtr_->dfdX(x, li, dfdx);
        }

        //- Return the jacobian w.r.t. the coefficients
        virtual void coeffJ
        (
            const typename UnivariateCoefficientEquation<Type>::XType& x,
            const label i,
            RectangularMatrix<Type>& J,
            const label li
        ) const
        {
            redirectEquationPtr_->coeffJ(x, i, J, li);
        }

        //- Set the pointer to the object registry
        virtual void setObr(const objectRegistry& obr)
        {
            UnivariateCoefficientEquation<Type>::setObr(obr);
            redirectEquationPtr_->setObr(obr);
        }

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const
        {
            writeCode(os);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CodedCoefficientEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
