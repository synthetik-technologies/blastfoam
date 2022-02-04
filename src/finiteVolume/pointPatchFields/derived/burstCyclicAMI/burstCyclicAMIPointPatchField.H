/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::burstCyclicAMIPointPatchField

Description
    Cyclic front and back plane patch field

SourceFiles
    burstCyclicAMIPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef burstCyclicAMIPointPatchField_H
#define burstCyclicAMIPointPatchField_H

#include "cyclicAMIPointPatchField.H"
#include "burstCyclicAMIPointPatch.H"
#include "burstPointPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class burstCyclicAMIPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class burstCyclicAMIPointPatchField
:
    public cyclicAMIPointPatchField<Type>,
    public burstPointPatchFieldBase
{
    // Private Data

        //- Intact point patch field
        autoPtr<pointPatchField<Type>> intactPointPatchField_;


public:

    //- Runtime type information
    TypeName(burstCyclicAMIPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        burstCyclicAMIPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        burstCyclicAMIPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        burstCyclicAMIPointPatchField
        (
            const burstCyclicAMIPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        burstCyclicAMIPointPatchField
        (
            const burstCyclicAMIPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new burstCyclicAMIPointPatchField<Type>
                (
                    *this, iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return the intact patch field
            const pointPatchField<Type>& intactPatchField() const
            {
                return intactPointPatchField_();
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const pointPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const pointPatchField<Type>&, const labelList&);


        // Evaluation functions

            //- Update coefficients
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        // Assignment operators

            virtual void operator=
            (
                const burstCyclicAMIPointPatchField<Type>&
            );

            virtual void operator=
            (
                const pointPatchField<Type>&
            );
            virtual void operator=(const Field<Type>&);
            virtual void operator=(const Type&);


        // Force an assignment irrespective of form of patch

            virtual void operator==
            (
                const burstCyclicAMIPointPatchField<Type>&
            );
            virtual void operator==
            (
                const pointPatchField<Type>&
            );
            virtual void operator==(const Field<Type>&);
            virtual void operator==(const Type&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "burstCyclicAMIPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
