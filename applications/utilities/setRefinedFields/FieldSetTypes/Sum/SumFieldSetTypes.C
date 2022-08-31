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

#include "SumFieldSetType.H"
#include "FieldSetTypesFwd.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldSetTypes
{

template<class Type>
using SumVol = Sum<Type, VolFieldSetType>;

template<class Type>
using SumSurface = Sum<Type, SurfaceFieldSetType>;

template<class Type>
using SumPoint = Sum<Type, PointFieldSetType>;
}

makeFieldSetTypeType(SumVol, scalar, VolFieldSetType);
makeFieldSetTypeType(SumVol, vector, VolFieldSetType);
makeFieldSetTypeType(SumVol, sphericalTensor, VolFieldSetType);
makeFieldSetTypeType(SumVol, symmTensor, VolFieldSetType);
makeFieldSetTypeType(SumVol, tensor, VolFieldSetType);

makeFieldSetTypeType(SumSurface, scalar, SurfaceFieldSetType);
makeFieldSetTypeType(SumSurface, vector, SurfaceFieldSetType);
makeFieldSetTypeType(SumSurface, sphericalTensor, SurfaceFieldSetType);
makeFieldSetTypeType(SumSurface, symmTensor, SurfaceFieldSetType);
makeFieldSetTypeType(SumSurface, tensor, SurfaceFieldSetType);

makeFieldSetTypeType(SumPoint, scalar, PointFieldSetType);
makeFieldSetTypeType(SumPoint, vector, PointFieldSetType);
makeFieldSetTypeType(SumPoint, sphericalTensor, PointFieldSetType);
makeFieldSetTypeType(SumPoint, symmTensor, PointFieldSetType);
makeFieldSetTypeType(SumPoint, tensor, PointFieldSetType);
}

// ************************************************************************* //

