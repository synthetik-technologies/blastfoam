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

\*---------------------------------------------------------------------------*/

#include "SumFieldSetType.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
Foam::FieldSetTypes::Sum<Type, FSType>::Sum
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const labelList& selectedIndices,
    Istream& is,
    const bool write
)
:
    FSType<Type>
    (
        mesh,
        dict,
        fieldName,
        selectedIndices,
        is,
        write
    )
{
    if (this->good_)
    {
        token t(is);
        if (t.isWord())
        {
            sumFld_.set(this->lookupOrRead(t.wordToken()));
        }
        else if
        (
            mesh.foundObject<FieldType>("sum:" + this->fieldPtr_->name())
        )
        {
            sumFld_.set
            (
                &mesh.lookupObject<FieldType>
                (
                    "sum:" + this->fieldPtr_->name()
                )
            );
        }
        else
        {
            is.putBack(t);
            Type sumVal = pTraits<Type>(is);
            FieldType* fldPtr =
                new FieldType
                (
                    IOobject
                    (
                        "sum:" + this->fieldPtr_->name(),
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimensioned<Type>
                    (
                        this->fieldPtr_->dimensions(),
                        sumVal
                    )
                );
            fldPtr->store(fldPtr);

            sumFld_.set
            (
                &mesh.lookupObject<FieldType>
                (
                    "sum:" + this->fieldPtr_->name()
                )
            );
        }

        wordList fldNames(is);
        flds_.resize(fldNames.size());
        forAll(flds_, i)
        {
            flds_.set (i, this->lookupOrRead(fldNames[i]));
        }

        this->setField();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
Foam::FieldSetTypes::Sum<Type, FSType>::~Sum()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
void Foam::FieldSetTypes::Sum<Type, FSType>::getInternalField
(
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{
    forAll(indices, i)
    {
        label celli = indices[i];
        Type res(sumFld_()[celli]);
        forAll(flds_, fldi)
        {
            res -= flds_[fldi][celli];
        }
        f[i] = res;
    }
}


template<class Type, template<class> class FSType>
void Foam::FieldSetTypes::Sum<Type, FSType>::getBoundaryField
(
    const label patchi,
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{
    forAll(indices, i)
    {
        label facei = indices[i];
        Type res(sumFld_().boundaryField()[patchi][facei]);
        forAll(flds_, fldi)
        {
            res -= flds_[fldi].boundaryField()[patchi][facei];
        }
        f[i] = res;
    }
}


// ************************************************************************* //
