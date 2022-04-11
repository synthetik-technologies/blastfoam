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
    Foam::TypeEquation

Description
    An equation using a scalarField as input, but outputing only a single,
    scalar value. Used in minimization.

\*---------------------------------------------------------------------------*/

#ifndef CodedUnivariateEquation_H
#define CodedUnivariateEquation_H

#include "univariateEquation.H"
#include "scalarField.H"
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
                          Class UnivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CodedUnivariateEquation
:
    public UnivariateEquation<Type>,
    public codedBase
{
    // Private Data

        //- Number of derivatives
        label nDerivatives_;

        //- The dynamically generated Equation pointer
        mutable autoPtr<univariateEquation<Type>> redirectEquationPtr_;


    // Private Member Functions

        //- Get the keywords associated with source code
        virtual wordList codeKeys() const;

        //- Adapt the context for the current object
        virtual void prepare
        (
            dynamicCode& dynCode,
            const dynamicCodeContext& context
        ) const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const;

        //- Compile, link and return the now coded Equation
        autoPtr<univariateEquation<Type>> compileNew();

        //- Expand strings in the code entry
        const dictionary& expandCodeDict(const dictionary& cDict) const;


public:

    TypeName("coded");

    //- Construct from dictionary
    CodedUnivariateEquation(const dictionary& dict);


    //- Destructor
    virtual ~CodedUnivariateEquation();


    // Member Functions

        //- Return the number of implemented derivatives
        virtual inline label nDerivatives() const
        {
            return redirectEquationPtr_->nDerivatives();
        }

        //- Return the function value
        virtual inline Type fX
        (
            const scalarList& x,
            const label li
        ) const
        {
            return redirectEquationPtr_->fX(x, li);
        }

        //- Return the function value
        virtual void dfdX
        (
            const scalarList& x,
            const label li,
            List<Type>& dfdx
        ) const
        {
            redirectEquationPtr_->dfdX(x, li, dfdx);
        }

        //- Return the function value
        virtual void d2fdX2
        (
            const scalarList& x,
            const label li,
            List<Type>& d2fdx2
        ) const
        {
            redirectEquationPtr_->d2fdX2(x, li, d2fdx2);
        }

        //- Return the function value
        virtual void d3fdX3
        (
            const scalarList& x,
            const label li,
            List<Type>& d3fdx3
        ) const
        {
            redirectEquationPtr_->d3fdX3(x, li, d3fdx3);
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
    #include "CodedUnivariateEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
