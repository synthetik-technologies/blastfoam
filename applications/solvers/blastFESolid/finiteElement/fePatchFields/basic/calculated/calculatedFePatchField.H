/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::calculatedFePatchField

Description
    A calculated boundary condition for feField

SourceFiles
    calculatedFePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef calculatedFePatchField_H
#define calculatedFePatchField_H

#include "fePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class calculatedFePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class calculatedFePatchField
:
    public fePatchField<Type>
{

public:

    //- Runtime type information
    TypeName("calculated");


    // Constructors

        //- Construct from patch and internal field
        calculatedFePatchField
        (
            const fePatch&,
            const DimensionedField<Type, feMesh>&
        );

        //- Construct from patch, internal field and dictionary
        calculatedFePatchField
        (
            const fePatch&,
            const DimensionedField<Type, feMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        calculatedFePatchField
        (
            const calculatedFePatchField<Type>&,
            const fePatch&,
            const DimensionedField<Type, feMesh>&,
            const fePatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        calculatedFePatchField
        (
            const calculatedFePatchField<Type>&,
            const DimensionedField<Type, feMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fePatchField<Type>>
        clone(const DimensionedField<Type, feMesh>& iF) const
        {
            return autoPtr<fePatchField<Type>>
            (
                new calculatedFePatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "calculatedFePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
