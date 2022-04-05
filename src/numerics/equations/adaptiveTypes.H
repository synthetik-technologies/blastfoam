/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#ifndef adaptiveTypes_H
#define adaptiveTypes_H

#include "scalar.H"
#include "vector.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class adaptiveError
{
public:

    template<class Type>
    static bool adaptive(const bool a)
    {
        return false;
    }

    template<class Type>
    static scalar error(const Type& f1, const Type& f2)
    {
        return great;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineAdaptiveType(Type)                                               \
    template<>                                                                 \
    bool adaptiveError::adaptive<Type>(const bool);                            \
                                                                               \
    template<>                                                                 \
    scalar adaptiveError::error<Type>(const Type&, const Type&);

#define defineAdaptiveIntegratorArgs(Dummy, Type)                              \
    template<>                                                                 \
    bool adaptiveError::adaptive<Type>(const bool);                            \
                                                                               \
    template<>                                                                 \
    scalar adaptiveError::error<Type>(const Type&, const Type&);

#define makeAdaptiveType(Func, Type)                                           \
    template<>                                                                 \
    bool adaptiveError::adaptive<Type>(const bool a) {return a;}               \
                                                                               \
    template<>                                                                 \
    scalar adaptiveError::error<Type>(const Type& f1, const Type& f2)          \
    {                                                                          \
        return Func;                                                           \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
