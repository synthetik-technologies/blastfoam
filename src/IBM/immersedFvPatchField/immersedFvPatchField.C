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

\*---------------------------------------------------------------------------*/

#include "immersedFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "immersedBoundaryObjectListSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedFvPatchField<Type>::immersedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF, Field<Type>(0))/*,
    object_
    (
        const_cast<immersedBoundaryObject&>
        (
            immersedBoundaryObjectListSolver::New
            (
                dynamicCast<const fvMesh>(p.boundaryMesh().mesh())
            ).objects()[p.patch().name()]
        )
    )*/
{}


template<class Type>
Foam::immersedFvPatchField<Type>::immersedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, Field<Type>(0))/*,
    object_
    (
        const_cast<immersedBoundaryObject&>
        (
            immersedBoundaryObjectListSolver::New
            (
                p.boundaryMesh().mesh()
            ).objects()[p.patch().name()]
        )
    )*/
{
    if (!isType<immersedFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::immersedFvPatchField<Type>::immersedFvPatchField
(
    const immersedFvPatchField<Type>&,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper&
)
:
    fvPatchField<Type>(p, iF, Field<Type>(0))/*,
    object_
    (
        const_cast<immersedBoundaryObject&>
        (
            immersedBoundaryObjectListSolver::New
            (
                p.boundaryMesh().mesh()
            ).objects()[p.patch().name()]
        )
    )*/
{
    if (!isType<immersedFvPatch>(p))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
Foam::immersedFvPatchField<Type>::immersedFvPatchField
(
    const immersedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf.patch(), iF, Field<Type>(0))/*,
    object_
    (
        const_cast<immersedBoundaryObject&>
        (
            immersedBoundaryObjectListSolver::New
            (
                ptf.patch().boundaryMesh().mesh()
            ).objects()[ptf.patch().patch().name()]
        )
    )*/
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// template<class Type>
// Foam::tmp<Foam::Field<Type>>
// Foam::immersedFvPatchField<Type>::patchInternalField() const
// {
//     immersedBoundaryObject& object
//     (
//         const_cast<immersedBoundaryObject&>
//         (
//             immersedBoundaryObjectListSolver::New
//             (
//                 dynamicCast<const fvMesh>(this->patch().boundaryMesh().mesh())
//             ).objects()[this->patch().patch().name()]
//         )
//     );
//     return object.patchExternalField(this->internalField());
// }
//
//
// template<class Type>
// void Foam::immersedFvPatchField<Type>::patchInternalField(Field<Type>& pif) const
// {
//     immersedBoundaryObject& object
//     (
//         const_cast<immersedBoundaryObject&>
//         (
//             immersedBoundaryObjectListSolver::New
//             (
//                 dynamicCast<const fvMesh>(this->patch().boundaryMesh().mesh())
//             ).objects()[this->patch().patch().name()]
//         )
//     );
//     pif = object.patchExternalField(this->internalField());
// }


template<class Type>
void Foam::immersedFvPatchField<Type>::updateCoeffs()
{
//     immersedBoundaryObject& object
//     (
//         const_cast<immersedBoundaryObject&>
//         (
//             immersedBoundaryObjectListSolver::New
//             (
//                 dynamicCast<const fvMesh>(this->patch().boundaryMesh().mesh())
//             ).objects()[this->patch().patch().name()]
//         )
//     );
//     object.setValues<Type>
//     (
//         const_cast<DimensionedField<Type, volMesh>&>
//         (
//             this->internalField()
//         )
//     );
}


// ************************************************************************* //
