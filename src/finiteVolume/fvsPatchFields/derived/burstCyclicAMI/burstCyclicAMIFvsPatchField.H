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
    Foam::burstCyclicAMIFvsPatchField

Description
    Foam::burstCyclicAMIFvsPatchField

SourceFiles
    burstCyclicAMIFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef burstCyclicAMIFvsPatchField_H
#define burstCyclicAMIFvsPatchField_H

#include "cyclicAMIFvsPatchField.H"
#include "burstCyclicAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class burstCyclicAMIFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class burstCyclicAMIFvsPatchField
:
    public cyclicAMIFvsPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(burstCyclicAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        burstCyclicAMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        burstCyclicAMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given burstCyclicAMIFvsPatchField onto a new patch
        burstCyclicAMIFvsPatchField
        (
            const burstCyclicAMIFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        burstCyclicAMIFvsPatchField(const burstCyclicAMIFvsPatchField<Type>&) = delete;

        //- Copy constructor setting internal field reference
        burstCyclicAMIFvsPatchField
        (
            const burstCyclicAMIFvsPatchField<Type>&,
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
                new burstCyclicAMIFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "burstCyclicAMIFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
