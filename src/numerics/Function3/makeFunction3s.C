/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
03-12-2021 Synthetik Applied Technologies : Added Function3
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

#include "None3.H"
#include "Constant3.H"
#include "ZeroConstant3.H"
#include "OneConstant3.H"
#include "Scale3.H"
#include "CodedFunction3.H"

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction3s(Type)                                                   \
    makeFunction3(Type);                                                       \
    makeFunction3Type(None, Type);                                             \
    makeFunction3Type(Constant, Type);                                         \
    makeFunction3Type(ZeroConstant, Type);                                     \
    makeFunction3Type(OneConstant, Type);                                      \
    makeFunction3Type(Scale, Type);                                            \
    makeFunction3Type(Coded, Type);

namespace Foam
{
    makeFunction3(label);
    makeFunction3Type(None, label);
    makeFunction3Type(Constant, label);

    makeFunction3s(scalar);
    makeFunction3s(vector);
    makeFunction3s(sphericalTensor);
    makeFunction3s(symmTensor);
    makeFunction3s(tensor);
}


// ************************************************************************* //
