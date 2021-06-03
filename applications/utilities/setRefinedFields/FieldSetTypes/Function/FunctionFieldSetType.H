/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::FieldSetTypes::Function

Description
    Linear combination of x, y, and z equations used by specifying an equation and the used coefficient functions i.e. f = a(x) + b(y)*d(x) + c(z)


Usage:
     volScalarFieldFunction p
    "-a(x) + b(y)*c(z)*d(x) - d(y)"
    {
        a
        {
            type uniform;
            value 1;
        }
        b
        {
            type uniform;
            value 0.5;
        }
        c
        {
            type uniform;
            value 1;
        }
        d
        {
            type uniform;
            value 0.5;
        }
    }


SourceFiles
    FunctionFieldSetType.C
    FunctionFieldSetType.C

\*---------------------------------------------------------------------------*/

#ifndef FunctionFieldSetType_H
#define FunctionFieldSetType_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldSetType.H"
#include "Function1.H"
#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Function Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class directionalFunction1
{
    //- Reference to Function1
    const Function1<Type>& func_;

    //- Component
    label cmpt_;


public:

    // Constructor
    directionalFunction1(const Function1<Type>& f, const label cmpt)
    :
        func_(f),
        cmpt_(cmpt)
    {}

    //- Destructor
    virtual ~directionalFunction1()
    {}


    // Public member functions

        //- Evaluate the function at a position
        Type value(const vector& xyz) const
        {
            return func_.value(xyz[cmpt_]);
        }

        //- Evaluate the function at a list of positions
        tmp<Field<Type>> value(const Field<vector>& xyz) const
        {
            tmp<Field<Type>> tmpf(new Field<Type>(xyz.size(), Zero));
            Field<Type>& f = tmpf.ref();
            forAll(f, i)
            {
                f[i] = func_.value(xyz[i][cmpt_]);
            }
            return tmpf;
        }
};


namespace FieldSetTypes
{
/*---------------------------------------------------------------------------*\
                          Class Function Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class Function
:
    public FieldSetType<Type, Patch, Mesh>
{
    //- Stored components of the equation
    //  Split by variables and operators
    tokenList equation_;

    //- List of Function1s, indexed by the coefficient name
    PtrListDictionary<Function1<Type>> Func1s_;

    //- List of functions and their corresponding components
    //  indexed by the coefficient name + the component, i.e. ax
    PtrListDictionary<directionalFunction1<Type>> funcs_;

    //- Set the selected cells
    virtual void setField();


public:

    TypeName("Function");

    // Constructor
    Function
    (
        const fvMesh& mesh,
        const word& fieldName,
        const labelList& selectedCells,
        Istream& is,
        const bool write
    );

    //- Destructor
    virtual ~Function();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FieldSetTypes
} // End namespace Foam

#ifdef NoRepository
    #include "FunctionFieldSetType.C"
#endif

#endif
// ************************************************************************* //
