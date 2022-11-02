/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
03-12-2021 Synthetik Applied Technologies : Added Function4
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "None4.H"
#include "Constant4.H"
#include "ZeroConstant4.H"
#include "OneConstant4.H"
#include "Scale4.H"
#include "CodedFunction4.H"

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction4s(Type)                                                   \
    makeFunction4(Type);                                                       \
    makeFunction4Type(None, Type);                                             \
    makeFunction4Type(Constant, Type);                                         \
    makeFunction4Type(ZeroConstant, Type);                                     \
    makeFunction4Type(OneConstant, Type);                                      \
    makeFunction4Type(Scale, Type);                                            \
    makeFunction4Type(Coded, Type);

namespace Foam
{
    makeFunction4(label);
    makeFunction4Type(None, label);
    makeFunction4Type(Constant, label);

    makeFunction4s(scalar);
    makeFunction4s(vector);
    makeFunction4s(sphericalTensor);
    makeFunction4s(symmTensor);
    makeFunction4s(tensor);
}


// ************************************************************************* //
