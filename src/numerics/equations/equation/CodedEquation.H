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

#ifndef CodedEquation_H
#define CodedEquation_H

#include "Equation.H"
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
class CodedEquation
:
    public Equation<Type>,
    public codedBase
{
    // Private Data

        //- Number of derivatives
        label nDerivatives_;

        //- The dynamically generated Equation pointer
        mutable autoPtr<equation<Type>> redirectEquationPtr_;


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
        autoPtr<equation<Type>> compileNew();

        //- Expand strings in the code entry
        const dictionary& expandCodeDict(const dictionary& cDict) const;


public:

    TypeName("coded");

    //- Construct from dictionary
    CodedEquation(const dictionary& dict);

    //- Destructor
    virtual ~CodedEquation();


    // Member Functions

        //- Return the number of implemented derivatives
        virtual inline label nDerivatives() const
        {
            return redirectEquationPtr_->nDerivatives();
        }

        //- Return the function value
        virtual inline Type fx
        (
            const scalar x,
            const label li
        ) const
        {
            return redirectEquationPtr_->fx(x, li);
        }

        //- Return the function value
        virtual inline Type dfdx
        (
            const scalar x,
            const label li
        ) const
        {
            return redirectEquationPtr_->dfdx(x, li);
        }

        //- Return the function value
        virtual inline Type d2fdx2
        (
            const scalar x,
            const label li
        ) const
        {
            return redirectEquationPtr_->d2fdx2(x, li);
        }

        //- Return the function value
        virtual inline Type d3fdx3
        (
            const scalar x,
            const label li
        ) const
        {
            return redirectEquationPtr_->d3fdx3(x, li);
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
    #include "CodedEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
