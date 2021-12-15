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
    Linear combination of x, y, and z equations used by specifying an equation
    and the used coefficient functions i.e. f = a(x) + b(y)*d(x) + c(z)


Usage:
     volScalarFieldFunction p
    {
        type constant;
        value 1.0
    }


SourceFiles
    FunctionFieldSetType.C
    FunctionFieldSetType.C

\*---------------------------------------------------------------------------*/

#ifndef FunctionFieldSetType_H
#define FunctionFieldSetType_H

#include "FieldSetType.H"
#include "Function3.H"
#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FieldSetTypes
{
/*---------------------------------------------------------------------------*\
                          Class Function Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class FSType>
class Function
:
    public FSType<Type>
{
    //- List of Function1s, indexed by the coefficient name
    autoPtr<Function3<Type>> func_;

    //- Return the subset of values given indices and positions
    virtual void getInternalField
    (
        const labelList& indices,
        const UIndirectList<vector>& pts,
        UIndirectList<Type>& f
    );

    //- Return the subset of values given indices and positions
    virtual void getBoundaryField
    (
        const label patchi,
        const labelList& indices,
        const UIndirectList<vector>& pts,
        UIndirectList<Type>& f
    );


public:

    TypeName("Function");

    // Constructor
    Function
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& fieldName,
        const labelList& selectedIndices,
        Istream& is,
        const bool write
    );

    //- Destructor
    virtual ~Function();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FieldSetTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FunctionFieldSetType.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
