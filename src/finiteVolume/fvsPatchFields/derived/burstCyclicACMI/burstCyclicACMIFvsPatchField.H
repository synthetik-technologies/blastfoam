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
    Foam::burstCyclicACMIFvsPatchField

Description
    Foam::burstCyclicACMIFvsPatchField

SourceFiles
    burstCyclicACMIFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef burstCyclicACMIFvsPatchField_H
#define burstCyclicACMIFvsPatchField_H

#include "cyclicACMIFvsPatchField.H"
#include "burstCyclicACMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class burstCyclicACMIFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class burstCyclicACMIFvsPatchField
:
    public cyclicACMIFvsPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(burstCyclicACMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        burstCyclicACMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        burstCyclicACMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given burstCyclicACMIFvsPatchField onto a new patch
        burstCyclicACMIFvsPatchField
        (
            const burstCyclicACMIFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        burstCyclicACMIFvsPatchField(const burstCyclicACMIFvsPatchField<Type>&) = delete;

        //- Copy constructor setting internal field reference
        burstCyclicACMIFvsPatchField
        (
            const burstCyclicACMIFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new burstCyclicACMIFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "burstCyclicACMIFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
