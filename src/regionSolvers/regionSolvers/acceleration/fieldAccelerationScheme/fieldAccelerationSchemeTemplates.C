/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "fieldAccelerationScheme.H"
#include "AccelerationSchemeBase.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template <class> class Patch, class Mesh>
Foam::fieldAccelerationScheme::fieldAccelerationScheme
(
    GeometricField<Type, Patch, Mesh>& field,
    const fvMesh& mesh,
    const labelList& patches,
    const dictionary& dict
)
:
    PtrList<accelerationScheme>(0),
    fieldName_(field.name()),
    mesh_(mesh)
{
    Info<< "Selecting acceleration schemes for " << field.name() << endl;
    incrIndent(Info);

    forAll(patches, pi)
    {
        const label patchi = patches[pi];
        PtrList<PatchFieldSelector<Type>> selectors
        (
            PatchFieldSelector<Type>::New
            (
                field.boundaryFieldRef()[patchi],
                dict.subOrEmptyDict(field.name())
            )
        );

        forAll(selectors, selectori)
        {
            const label acceleratori = this->size();
            this->setSize(acceleratori + 1);
            this->set
            (
                acceleratori,
                AccelerationSchemeBase<Type, Patch, Mesh>::New
                (
                    field,
                    selectors.set(selectori, nullptr),
                    dict
                )
            );
        }
    }
    this->read(dict);

    Info<< decrIndent << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template <class> class Patch, class Mesh>
bool Foam::fieldAccelerationScheme::setPatchAccelerators
(
    const labelList& patches,
    const dictionary& dict
)
{
    if (mesh_.foundObject<GeometricField<Type, Patch, Mesh>>(fieldName_))
    {
        GeometricField<Type, Patch, Mesh>& field =
            mesh_.lookupObjectRef<GeometricField<Type, Patch, Mesh>>
            (
                fieldName_
            );
        forAll(patches, pi)
        {
            const label patchi = patches[pi];
            PtrList<PatchFieldSelector<Type>> selectors
            (
                PatchFieldSelector<Type>::New
                (
                    field.boundaryFieldRef()[patchi],
                    dict.subOrEmptyDict(fieldName_)
                )
            );

            forAll(selectors, selectori)
            {
                const label acceleratori = this->size();
                this->setSize(acceleratori + 1);
                this->set
                (
                    acceleratori,
                    AccelerationSchemeBase<Type, Patch, Mesh>::New
                    (
                        field,
                        selectors.set(selectori, nullptr),
                        dict
                    )
                );
            }
        }
        return true;
    }
    return false;
}


template<class Type, template <class> class Patch, class Mesh>
bool Foam::fieldAccelerationScheme::storePrevIter()
{
    if (mesh_.foundObject<GeometricField<Type, Patch, Mesh>>(fieldName_))
    {
        //- Store the previous iterator for the field
        mesh_.lookupObjectRef<GeometricField<Type, Patch, Mesh>>
        (
            fieldName_
        ).storePrevIter();


        // Set all patches as needed (i.e. gradient)
        forAll(*this, i)
        {
            this->operator[](i).storePrevIter();
        }
        return true;
    }
    return false;
}

// ************************************************************************* //
