/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "volFields.H"
#include "cellCellStencil.H"
#include "cellCellStencilObject.H"
#include "dynamicOversetFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(p, iF),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchField<Type>(ptf, p, iF, mapper),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchField<Type>(p, iF, dict),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(ptf, iF),
    oversetPatch_(ptf.oversetPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::oversetFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (this->oversetPatch_.master())
    {
        // Trigger interpolation
        const fvMesh& mesh = this->internalField().mesh();
        const dictionary& fvSchemes = mesh.schemesDict();
        const word& fldName = this->internalField().name();

        if (&mesh.lduAddr() != &mesh.fvMesh::lduAddr())
        {
            // Running extended addressing. Flux correction already done
            // in the linear solver (through the initUpdateInterfaceMatrix
            // method below)
            if (debug)
            {
                Info<< "Skipping overset interpolation for solved-for field "
                    << fldName << endl;
            }
        }
        else if (!fvSchemes.found("oversetInterpolation"))
        {
            IOWarningInFunction(fvSchemes)
                << "Missing required dictionary entry"
                << " 'oversetInterpolation'"
                << ". Skipping overset interpolation for field "
                << fldName << endl;
        }
        else if (fvSchemes.found("oversetInterpolationRequired"))
        {
            // Backwards compatibility mode: only interpolate what is
            // explicitly mentioned

            if (fvSchemes.found("oversetInterpolationSuppressed"))
            {
                FatalIOErrorInFunction(fvSchemes)
                    << "Cannot have both dictionary entry"
                    << " 'oversetInterpolationSuppresed' and "
                    << " 'oversetInterpolationRequired' for field "
                    << fldName << exit(FatalIOError);
            }
            const dictionary& intDict = fvSchemes.subDict
            (
                "oversetInterpolationRequired"
            );
            if (intDict.found(fldName))
            {
                if (debug)
                {
                    Info<< "Interpolating field " << fldName << endl;
                }

                // Interpolate without bc update (since would trigger infinite
                // recursion back to oversetFvPatchField<Type>::evaluate)
                // The only problem is bcs that use the new cell values
                // (e.g. zeroGradient, processor). These need to appear -after-
                // the 'overset' bc.
                dynamicCast<const dynamicOversetFvMesh>(mesh).interpolate
                (
                    const_cast<Field<Type>&>
                    (
                        this->primitiveField()
                    )
                );
            }
            else if (debug)
            {
                Info<< "Skipping overset interpolation for field "
                    << fldName << endl;
            }
        }
        else
        {
            const dictionary* dictPtr
            (
                fvSchemes.subDictPtr("oversetInterpolationSuppressed")
            );

            const wordHashSet& suppress =
                Stencil::New(mesh).nonInterpolatedFields();

            bool skipInterpolate = suppress.found(fldName);

            if (dictPtr)
            {
                skipInterpolate =
                    skipInterpolate
                 || dictPtr->found(fldName);
            }

            if (skipInterpolate)
            {
                if (debug)
                {
                    Info<< "Skipping suppressed overset interpolation"
                        << " for field " << fldName << endl;
                }
            }
            else
            {
                if (debug)
                {
                    Info<< "Interpolating non-suppressed field " << fldName
                        << endl;
                }
                dynamicCast<const dynamicOversetFvMesh>(mesh).interpolate
                (
                    const_cast<Field<Type>&>
                    (
                        this->primitiveField()
                    )
                );
            }
        }
    }
}


template<class Type>
void Foam::oversetFvPatchField<Type>::write(Ostream& os) const
{
    zeroGradientFvPatchField<Type>::write(os);
    // Make sure to write the value for ease of postprocessing.
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
