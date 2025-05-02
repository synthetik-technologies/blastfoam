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

#include "fePatchField.H"
#include "feMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fePatchField<Type>::fePatchField
(
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF
)
:
    patch_(p),
    internalField_(iF),
    updated_(false)
{}


template<class Type>
Foam::fePatchField<Type>::fePatchField
(
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF,
    const dictionary& dict
)
:
    patch_(p),
    internalField_(iF),
    updated_(false)
{}


template<class Type>
Foam::fePatchField<Type>::fePatchField
(
    const fePatchField<Type>& ptf,
    const fePatch& p,
    const DimensionedField<Type, feMesh>& iF,
    const fePatchFieldMapper&
)
:
    patch_(p),
    internalField_(iF),
    updated_(false)
{}


template<class Type>
Foam::fePatchField<Type>::fePatchField
(
    const fePatchField<Type>& ptf,
    const DimensionedField<Type, feMesh>& iF
)
:
    patch_(ptf.patch_),
    internalField_(iF),
    updated_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::objectRegistry& Foam::fePatchField<Type>::db() const
{
    return patch_.boundaryMesh().mesh().thisDb();
}


template<class Type>
void Foam::fePatchField<Type>::write(Ostream& os) const
{
    writeEntry(os, "type", type());

    if (overridesConstraint())
    {
        writeEntry(os, "patchType", patch().type());
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fePatchField<Type>::patchInternalField() const
{
    return patchInternalField(primitiveField());
}


template<class Type>
template<class Type1>
Foam::tmp<Foam::Field<Type1>>
Foam::fePatchField<Type>::patchInternalField
(
    const Field<Type1>& iF,
    const labelList& meshFes
) const
{
    // Check size
    if (iF.size() != primitiveField().size())
    {
        FatalErrorInFunction
            << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << primitiveField().size()
            << abort(FatalError);
    }

    return tmp<Field<Type1>>(new Field<Type1>(iF, meshFes));
}


template<class Type>
template<class Type1>
Foam::tmp<Foam::Field<Type1>>
Foam::fePatchField<Type>::patchInternalField
(
    const Field<Type1>& iF
) const
{
    return patchInternalField(iF, patch().meshNodes());
}


template<class Type>
template<class Type1>
void Foam::fePatchField<Type>::addToInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF
) const
{
    // Check size
    if (iF.size() != primitiveField().size())
    {
        FatalErrorInFunction
            << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << primitiveField().size()
            << abort(FatalError);
    }

    if (pF.size() != size())
    {
        FatalErrorInFunction
            << "given patch field does not correspond to the mesh. "
            << "Field size: " << pF.size()
            << " mesh size: " << size()
            << abort(FatalError);
    }

    // Get the addressing
    const labelList& mp = patch().meshNodes();

    forAll(mp, fei)
    {
        iF[mp[fei]] += pF[fei];
    }
}


template<class Type>
template<class Type1>
void Foam::fePatchField<Type>::addToInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF,
    const labelList& fes
) const
{
    // Check size
    if (iF.size() != primitiveField().size())
    {
        FatalErrorInFunction
            << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << primitiveField().size()
            << abort(FatalError);
    }

    if (pF.size() != size())
    {
        FatalErrorInFunction
            << "given patch field does not correspond to the mesh. "
            << "Field size: " << pF.size()
            << " mesh size: " << size()
            << abort(FatalError);
    }

    // Get the addressing
    const labelList& mp = patch().meshNodes();

    forAll(fes, i)
    {
        label fei = fes[i];
        iF[mp[fei]] += pF[fei];
    }
}


template<class Type>
template<class Type1>
void Foam::fePatchField<Type>::setInInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF,
    const labelList& meshNodes
) const
{
    // Check size
    if (iF.size() != primitiveField().size())
    {
        FatalErrorInFunction
            << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << primitiveField().size()
            << abort(FatalError);
    }

    if (pF.size() != meshNodes.size())
    {
        FatalErrorInFunction
            << "given patch field does not correspond to the meshNodes. "
            << "Field size: " << pF.size()
            << " meshNodes size: " << size()
            << abort(FatalError);
    }

    forAll(meshNodes, fei)
    {
        iF[meshNodes[fei]] = pF[fei];
    }
}


template<class Type>
template<class Type1>
void Foam::fePatchField<Type>::setInInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF
) const
{
    setInInternalField(iF, pF, patch().meshNodes());
}


template<class Type>
void Foam::fePatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated_)
    {
        updateCoeffs();
    }

    updated_ = false;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const fePatchField<Type>& ptf
)
{
    ptf.write(os);

    os.check("Ostream& operator<<(Ostream&, const fePatchField<Type>&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fePatchFieldNew.C"

// ************************************************************************* //
