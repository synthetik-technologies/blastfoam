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

SourceFiles


\*---------------------------------------------------------------------------*/

#include "InitialValueFieldSetType.H"
#include "FieldSetTypesFwd.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldSetTypes
{

template<class Type>
using InitialValueVol = InitialValue<Type, VolFieldSetType>;

template<class Type>
using InitialValueSurface = InitialValue<Type, SurfaceFieldSetType>;

template<class Type>
using InitialValuePoint = InitialValue<Type, PointFieldSetType>;
}

makeFieldSetTypeType(InitialValueVol, scalar, VolFieldSetType);
makeFieldSetTypeType(InitialValueVol, vector, VolFieldSetType);
makeFieldSetTypeType(InitialValueVol, sphericalTensor, VolFieldSetType);
makeFieldSetTypeType(InitialValueVol, symmTensor, VolFieldSetType);
makeFieldSetTypeType(InitialValueVol, tensor, VolFieldSetType);

makeFieldSetTypeType(InitialValueSurface, scalar, SurfaceFieldSetType);
makeFieldSetTypeType(InitialValueSurface, vector, SurfaceFieldSetType);
makeFieldSetTypeType(InitialValueSurface, sphericalTensor, SurfaceFieldSetType);
makeFieldSetTypeType(InitialValueSurface, symmTensor, SurfaceFieldSetType);
makeFieldSetTypeType(InitialValueSurface, tensor, SurfaceFieldSetType);

makeFieldSetTypeType(InitialValuePoint, scalar, PointFieldSetType);
makeFieldSetTypeType(InitialValuePoint, vector, PointFieldSetType);
makeFieldSetTypeType(InitialValuePoint, sphericalTensor, PointFieldSetType);
makeFieldSetTypeType(InitialValuePoint, symmTensor, PointFieldSetType);
makeFieldSetTypeType(InitialValuePoint, tensor, PointFieldSetType);
}
// ************************************************************************* //

