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

#include "fieldAccelerationSchemeList.H"
#include "globalPolyBoundaryMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::fieldAccelerationSchemeList::addField
(
    GeometricField<Type, Patch, Mesh>& field
)
{
    addField(field, globalPolyBoundaryMesh::New(mesh_).coupledPatches());
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::fieldAccelerationSchemeList::addField
(
    GeometricField<Type, Patch, Mesh>& field,
    const labelList& patches
)
{
    if (this->found(field.name()) || !dict_.isDict(subDictName_))
    {
        return;
    }

    const label fieldi = this->size();
    this->setSize(fieldi + 1);
    this->set
    (
        fieldi,
        field.name(),
        new fieldAccelerationScheme
        (
            field,
            mesh_,
            patches,
            (subDictName_ != word::null)
          ? dict_.subDict(subDictName_)
          : dict_
        )
    );
}


// ************************************************************************* //
