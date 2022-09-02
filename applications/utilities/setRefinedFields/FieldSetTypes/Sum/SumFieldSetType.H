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
    Foam::FieldSetTypes::Sum

Description


Usage:
    Set field values based on the sum of the list of fields provided and the
    total value/fiel provided

    fieldValues
    (
        volScalarFieldSum alpha.water 1.0 (alpha.air alpha.solid)
        volScalarFieldSum Y1.water sumY (Y.air Y.solid)
    );

SourceFiles
    SumFieldSetType.C
    SumFieldSetType.C

\*---------------------------------------------------------------------------*/

#ifndef SumFieldSetType_H
#define SumFieldSetType_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldSetType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldSetTypes
{

/*---------------------------------------------------------------------------*\
                          Class Sum Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class FSType>
class Sum
:
    public FSType<Type>
{
    typedef typename FSType<Type>::FieldType FieldType;

    //- Sum field
    UautoPtr<const FieldType> sumFld_;

    //- Uniform value of sum of fields
    Type sum_;

    //- List of fields to sum
    UPtrList<const FieldType> flds_;

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

    TypeName("Sum");

    // Constructor
    Sum
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& fieldName,
        const labelList& selectedIndices,
        Istream& is,
        const bool write
    );

    //- Destructor
    virtual ~Sum();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FieldSetTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SumFieldSetType.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
