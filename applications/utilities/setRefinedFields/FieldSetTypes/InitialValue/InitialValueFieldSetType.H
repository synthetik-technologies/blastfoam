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
    Foam::FieldSetTypes::InitialValue

Description
    Use the saved field value as the field value. This is used to reset a field
    to the initial conditions

SourceFiles
    InitialValueFieldSetType.C
    InitialValueFieldSetType.C

\*---------------------------------------------------------------------------*/

#ifndef InitialValueFieldSetType_H
#define InitialValueFieldSetType_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldSetType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldSetTypes
{
/*---------------------------------------------------------------------------*\
                          Class InitialValue Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class InitialValue
:
    public FieldSetType<Type, Patch, Mesh>
{
    //- Value to set
    UautoPtr<GeometricField<Type, Patch, Mesh>> origFieldPtr_;

    //- Set the selected cells
    virtual void setField();


public:

    TypeName("InitialValue");

    // Constructor
    InitialValue
    (
        const fvMesh& mesh,
        const word& fieldName,
        const labelList& selectedCells,
        Istream& is,
        const bool write
    );

    //- Destructor
    virtual ~InitialValue();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FieldSetTypes
} // End namespace Foam

#ifdef NoRepository
    #include "InitialValueFieldSetType.C"
#endif

#endif
// ************************************************************************* //
