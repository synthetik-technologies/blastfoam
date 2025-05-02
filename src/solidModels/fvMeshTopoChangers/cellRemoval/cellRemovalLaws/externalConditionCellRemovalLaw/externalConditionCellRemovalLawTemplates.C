/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "externalConditionCellRemovalLaw.H"


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::externalConditionCellRemovalLaw::saveGeoFields
(
    objectRegistry& obr
) const
{
    typedef GeometricField<Type, Patch, Mesh> GeoField;
    HashTable<GeoField*> fields = obr.lookupClass<GeoField>();

    forAllIter
    (
        typename HashTable<GeoField*>,
        fields,
        iter
    )
    {
        const GeoField& field = *iter();
        GeoField* interpField = subsetter_->interpolate(*iter()).ptr();
        interpField->rename(field.name());
        interpField->writeOpt() = field.writeOpt();
        interpField->store(interpField);
    }
}


template<class Type, class CombineOp>
void Foam::externalConditionCellRemovalLaw::map
(
    const word& fieldName,
    const CombineOp& cop,
    DimensionedField<Type, volMesh>& field
) const
{
    if (!active())
    {
        return;
    }
    const GeometricField<Type, fvPatchField, volMesh>& rfield =
        subsetter_->subMesh().lookupObject
        <
            GeometricField<Type, fvPatchField, volMesh>
        >(fieldName);
    meshMap().mapTgtToSrc(rfield, cop, field);
}

// ************************************************************************* //
