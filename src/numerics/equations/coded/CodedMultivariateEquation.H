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
    Foam::CodedMultvariateEquation

Description
    User defined equation with multiple inputs and multiple output

\*---------------------------------------------------------------------------*/

#ifndef CodedMultvariateEquation_H
#define CodedMultvariateEquation_H

#include "regEquation.H"
#include "MultivariateEquation.H"
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
                    Class CodedMultvariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CodedMultivariateEquation
:
    public regEquation<Type, MultivariateEquation>,
    public codedBase
{
    // Private Data

        //- Number of derivatives
        label nDerivatives_;

        //- The dynamically generated Equation pointer
        mutable autoPtr<multivariateEquation<Type>> redirectEquationPtr_;


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
        autoPtr<multivariateEquation<Type>> compileNew();

        //- Expand strings in the code entry
        const dictionary& expandCodeDict(const dictionary& cDict) const;


public:

    // Constructors

        //- Construct from name and dictionary
        CodedMultivariateEquation
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

    //- Destructor
    virtual ~CodedMultivariateEquation();


    // Member Functions

        //- Return the number of implemented derivatives
        virtual inline label nDerivatives() const
        {
            return redirectEquationPtr_->nDerivatives();
        }

        //- Return the function value
        virtual inline void FX
        (
            const UList<scalar>& x,
            const label li,
            List<Type>& fx
        ) const
        {
            redirectEquationPtr_->FX(x, li, fx);
        }

        //- Return the function value
        virtual void jacobian
        (
            const UList<scalar>& x,
            const label li,
            List<Type>& fx,
            RectangularMatrix<Type>& J
        ) const
        {
            redirectEquationPtr_->jacobian(x, li, fx, J);
        }

        //- Write data to dictionary stream
        using regIOobject::write;
        virtual void write(Ostream& os) const
        {
            writeCode(os);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CodedMultivariateEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
