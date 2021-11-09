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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::immersedFvPatchField<Type>::setInternalField
(
    const Type& val
)
{

    object_.setInternal
    (
        this->internalFieldRef(),
        val
    );
}


template<class Type>
void Foam::immersedFvPatchField<Type>::setInternalField
(
    const Field<Type>& vals
)
{

    object_.setInternal
    (
        this->internalFieldRef(),
        vals
    );
}


template<class Type>
void Foam::immersedFvPatchField<Type>::setInternalFieldAvg()
{
    Field<Type>& iF(this->internalFieldRef());
    const labelList& internalC(object_.internalCells());
    const scalarField& V = iF.mesh().V();
    Type sumValV(Zero);
    scalar sumV(0.0);

    forAll(internalC, i)
    {
        const label celli = internalC[i];
        sumV += V[celli];
        sumValV += V[celli]*iF[celli];
    }

    object_.setInternal
    (
        iF,
        returnReduce(sumValV, sumOp<Type>())
       /returnReduce(sumV, sumOp<scalar>())
    );
}


template<class Type>
void Foam::immersedFvPatchField<Type>::setPatchFieldAvg
(
    const Field<Type>& pf
)
{
    const scalarField& magSf = object_.magSf();
    const labelList& pI = object_.patchInternalCells();
    Type sumValSf(Zero);
    scalar sumSf(0.0);

    forAll(pI, i)
    {
        label celli = pI[i];
        if (celli >= 0)
        {
            sumValSf += magSf[i];
            sumSf += magSf[i]*pf[i];
        }
    }

    if (sumSf < small)
    {
        internalFieldRef() = Zero;
        return;
    }

    object_.setInternal
    (
        internalFieldRef(),
        returnReduce(sumValSf, sumOp<Type>())
       /returnReduce(sumSf, sumOp<Type>())
    );
}


template<class Type>
void Foam::immersedFvPatchField<Type>::setPatchInternalField
(
    const Field<Type>& pf
)
{
    const scalarField& magSf = object_.magSf();
    const labelList& pI = object_.patchInternalCells();
    Field<Type> iFSf(object_.size(), Zero);
    scalarField sumMagSf(object_.size(), Zero);
    globalIndex globalCells
    (
        this->patch().boundaryMesh().mesh().nCells()
    );

    label npI = 0;
    Map<label> map;
    forAll(pI, i)
    {
        label celli = pI[i];
        if (celli >= 0)
        {
            label gCelli = globalCells.toGlobal(celli);
            if (!map.found(gCelli))
            {
                map.insert(gCelli, npI++);
            }
            label fi = map[gCelli];

            iFSf[fi] += magSf[i]*pf[i];
            sumMagSf[fi] += magSf[i];
        }
    }
    reduce(iFSf, sumOp<List<Type>>());
    reduce(sumMagSf, sumOp<List<scalar>>());
    iFSf /= max(sumMagSf, small);

    Field<Type>& iF(internalFieldRef());
    forAllConstIter(Map<label>, map, iter)
    {
        label gCelli = iter.key();
        if (globalCells.isLocal(gCelli))
        {
            label fi = iter();
            iF[globalCells.toLocal(gCelli)] = iFSf[fi];
        }
    }
}


template<class Type>
void Foam::immersedFvPatchField<Type>::smoothInternalField
(
    Field<Type>& iF
)
{
    const labelList& internalC = this->object_.internalCells();
    const fvMesh& mesh = this->internalField().mesh();
    const extendedNLevelCPCCellToCellStencil& stencil =
        extendedNLevelCPCCellToCellStencil::New(mesh);
    Field<Type> newValues(internalC.size(), Zero);
    const scalarField& V = mesh.V();
    List<List<Type>> stencilFld(mesh.nCells());
    List<List<scalar>> VFld(mesh.nCells());
    stencil.collectData(iF, stencilFld);
    stencil.collectData(V, VFld);

    forAll(internalC, i)
    {
        label celli = internalC[i];
        scalar sumV = 0.0;
        forAll(VFld[celli], j)
        {
            sumV += VFld[celli][j];
            newValues[i] += VFld[celli][j]*stencilFld[celli][j];
        }
        newValues[i] /= sumV;
    }

    forAll(internalC, i)
    {
        iF[internalC[i]] = newValues[i];
    }
}


template<class Type>
Foam::Field<Type>&
Foam::immersedFvPatchField<Type>::internalFieldRef()
{
    const Field<Type>& iF(this->internalField());
    return const_cast<Field<Type>&>(iF);
}


template<class Type>
void Foam::immersedFvPatchField<Type>::setInternalCells
(
    fvMatrix<Type>& m
) const
{
    // Fix the value in dead cells
    if (setInternal_)
    {
        const labelList& cells = object_.allInternalCells();

        // Boost the diagonal of dead cells by the volume ratio
        // Volume ratio is set to SMALL; revert for diagonal
        // This should also guarantee strong diagonal dominance.
        // HJ, 19/Jun/2019

        scalarField& diag = m.diag();

        forAll (cells, celli)
        {
            diag[cells[celli]] = great;
        }

        // Set values
        m.setValues
        (
            cells,
            Field<Type>(cells.size(), internalValue_)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedFvPatchField<Type>::immersedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF, Field<Type>(0)),
    object_
    (
        dynamicCast<const immersedPolyPatch&>
        (
            p.patch()
        ).immersedObject()
    ),
    setPatchInternal_(false),
    setInternal_(false),
    internalValue_(Zero),
    nSmooth_(0)
{}


template<class Type>
Foam::immersedFvPatchField<Type>::immersedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, Field<Type>(0)),
    object_
    (
        dynamicCast<const immersedPolyPatch&>
        (
            p.patch()
        ).immersedObject()
    ),
    setPatchInternal_(dict.lookupOrDefault("setPatchInternal", false)),
    setInternal_(dict.lookupOrDefault("setInternal", false)),
    internalValue_(Zero),
    nSmooth_(dict.lookupOrDefault<label>("nSmooth", 0))
{
    if (!isType<immersedFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << p.type() << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::immersedFvPatchField<Type>::immersedFvPatchField
(
    const immersedFvPatchField<Type>& ip,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper&
)
:
    fvPatchField<Type>(p, iF, Field<Type>(0)),
    object_
    (
        dynamicCast<const immersedPolyPatch&>
        (
            p.patch()
        ).immersedObject()
    ),
    setPatchInternal_(ip.setPatchInternal_),
    setInternal_(ip.setInternal_),
    internalValue_(ip.internalValue_),
    nSmooth_(ip.nSmooth_)
{
    if (!isType<immersedFvPatch>(p))
    {
        FatalErrorInFunction
            << "' not constraint type '" << p.type() << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
Foam::immersedFvPatchField<Type>::immersedFvPatchField
(
    const immersedFvPatchField<Type>& ip,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ip.patch(), iF, Field<Type>(0)),
    object_
    (
        dynamicCast<const immersedPolyPatch&>
        (
            ip.patch().patch()
        ).immersedObject()
    ),
    setPatchInternal_(ip.setPatchInternal_),
    setInternal_(ip.setInternal_),
    internalValue_(ip.internalValue_),
    nSmooth_(ip.nSmooth_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type>
// Foam::tmp<Foam::Field<Type>>
// Foam::immersedFvPatchField<Type>::patchInternalField() const
// {
//     return object_.patchExternalField(this->internalField());
// }
//
//
// template<class Type>
// void Foam::immersedFvPatchField<Type>::patchInternalField(Field<Type>& pif) const
// {
//     pif = object_.patchExternalField(this->internalField());
// }


template<class Type>
void Foam::immersedFvPatchField<Type>::updateCoeffs()
{
    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
void
Foam::immersedFvPatchField<Type>::addForcing
(
    Field<Type>& F,
    const Field<scalar>& alphaRho,
    const Field<Type>& old,
    const Field<Type>& RHS,
    const scalar& dt
) const
{}


template<class Type>
void Foam::immersedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (setPatchInternal_)
    {
        writeEntry(os, "setPatchInternal", setPatchInternal_);
    }
    if (setInternal_)
    {
        writeEntry(os, "setInternal", setInternal_);
        writeEntry(os, "internalValue", internalValue_);
    }
}


// ************************************************************************* //
