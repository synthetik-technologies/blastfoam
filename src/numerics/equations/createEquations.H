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

#ifndef createEquations_H
#define createEquations_H
#include "equation.H"
#include "UnivariateEquation.H"

// * * * * * * * * * * * * * * Helper MACROS * * * * * * * * * * * * * * * * //

#define VA_NARGS_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define VA_NARGS(...) VA_NARGS_IMPL(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)

// Make a FOREACH macro
#define FE_0(F, X, Y) F(X,Y)
#define FE_1(F, X, Y, ...) F(X,Y) FE_0(F, X, __VA_ARGS__)
#define FE_2(F, X, Y, ...) F(X,Y) FE_1(F, X, __VA_ARGS__)
#define FE_3(F, X, Y, ...) F(X,Y) FE_2(F, X, __VA_ARGS__)
//... repeat as needed

#define GET_MACRO(_0,_1,_2,_3,NAME,...) NAME
#define FOR_EACH(F,X,...) \
  GET_MACRO(__VA_ARGS__,FE_3,FE_2,FE_1,FE_0)(F,X,__VA_ARGS__)


// * * * * * * * * * * * * * * Define Functions  * * * * * * * * * * * * * * //

#define defineEquation(Type, FuncName)                                         \
    virtual Type FuncName(const scalar, const label) const;                    \
    virtual string print##FuncName() const;

#define defineEquationClass(Name, Base, Type, Min, Max, ...)                   \
class Name                                                                     \
:                                                                              \
    public Base                                                                \
{                                                                              \
public:                                                                        \
    Name() : Base(Min, Max) {}                                                 \
                                                                               \
    virtual ~Name() {}                                                         \
                                                                               \
    virtual label nDerivatives() const                                         \
    {                                                                          \
        return VA_NARGS(__VA_ARGS__);                                          \
    }                                                                          \
    FOR_EACH(defineEquation, Type, __VA_ARGS__);                              \
};

#define defineEquation0(Name, Min, Max)                                        \
defineEquationClass(Name, equation, scalar, Min, Max, fx)

#define defineEquation1(Name, Min, Max)                                        \
defineEquationClass(Name, equation, scalar, Min, Max, fx, dfdx)

#define defineEquation2(Name, Min, Max)                                        \
defineEquationClass(Name, equation, scalar, Min, Max, fx, dfdx, d2fdx2)

#define defineEquation3(Name, Min, Max)                                        \
defineEquationClass(Name, equation, scalar, Min, Max, fx, dfdx, d2fdx2, d3fdx3)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineUnivariateEquation0(Name, Type, Min, Max)                        \
defineEquationClass(Name, UnivariateEquation<Type>, Type, Min, Max, fx)

#define defineUnivariateEquation1(Name, Type, Min, Max)                        \
defineEquationClass(Name, UnivariateEquation<Type>, Type, Min, Max, fx, dfdx)

#define defineUnivariateEquation2(Name, Type, Min, Max)                        \
defineEquationClass(Name, UnivariateEquation<Type>, Type, Min, Max, fx, dfdx, d2fdx2)

#define defineUnivariateEquation3(Name, Type, Min, Max)                        \
defineEquationClass(Name, UnivariateEquation<Type>, Type, Min, Max, fx, dfdx, d2fdx2, d3fdx3)



// * * * * * * * * * * * * * * * Add Functions * * * * * * * * * * * * * * * //

#define addEquationFunc(Name, Type, FuncName, Func)                            \
Type Name::FuncName(const scalar x, const label li) const {return Func;}       \
string Name::print##FuncName() const {return #FuncName "(x) = " #Func;}

#define addUnivariateEquation0Funcs(Name, Type, FX)                            \
addEquationFunc(Name, Type, fx, FX)

#define addUnivariateEquation1Funcs(Name, Type, FX, dFdX)                      \
addEquationFunc(Name, Type, fx, FX)                                            \
addEquationFunc(Name, Type, dfdx, dFdX)

#define addUnivariateEquation2Funcs(Name, Type, FX, dFdX, d2FdX2)              \
addEquationFunc(Name, Type, fx, FX)                                            \
addEquationFunc(Name, Type, dfdx, dFdX)                                        \
addEquationFunc(Name, Type, d2fdx2, d2FdX2)

#define addUnivariateEquation3Funcs(Name, Type, FX, dFdX, d2FdX2, d3FdX3)      \
addEquationFunc(Name, Type, fx, FX)                                            \
addEquationFunc(Name, Type, dfdx, dFdX)                                        \
addEquationFunc(Name, Type, d2fdx2, d2FdX2)                                    \
addEquationFunc(Name, Type, d3fdx3, d3FdX3)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addEquation0Funcs(Name, FX)                                            \
addUnivariateEquation0Funcs(Name, scalar, FX)

#define addEquation1Funcs(Name, FX, dFdX)                                      \
addUnivariateEquation1Funcs(Name, scalar, FX, dFdX)

#define addEquation2Funcs(Name, FX, dFdX, d2FdX2)                              \
addUnivariateEquation2Funcs(Name, scalar, FX, dFdX, d2FdX2)

#define addEquation3Funcs(Name, FX, dFdX, d2FdX2, d3FdX3)                      \
addUnivariateEquation3Funcs(Name, scalar, FX, dFdX, d2FdX2, d3DdX3)


// * * * * * * * * * * * * Define and Implement Functions  * * * * * * * * * //

#define createEquation0(Name, Min, Max, FX)                                    \
defineEquation0(Name, Min, Max)                                                \
addEquation0Funcs(Name, FX)

#define createEquation1(Name, Min, Max, FX, dFdX)                              \
defineEquation1(Name, Min, Max)                                                \
addEquation1Funcs(Name, FX, dFdX)

#define createEquation2(Name, Min, Max, FX, dFdX, d2FdX2)                      \
defineEquation2(Name, Min, Max)                                                \
addEquation2Funcs(Name, FX, dFdX, d2FdX2)

#define createEquation3(Name, Min, Max, FX, dFdX, d2FdX2, d3FdX3)              \
defineEquation3(Name, Min, Max)                                                \
addEquation3Funcs(Name, FX, dFdX, d2FdX2, d3FdX3)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define createUnivariateEquation0(Name, Type, Min, Max, FX)                    \
defineUnivariateEquation0(Name, Type, Min, Max)                                \
addUnivariateEquation0Funcs(Name, Type, FX)

#define createUnivariateEquation1(Name, Type, Min, Max, FX, dFdX)              \
defineUnivariateEquation1(Name, Type, Min, Max)                                \
addUnivariateEquation1Funcs(Name, Type, FX, dFdX)

#define createUnivariateEquation2(Name, Type, Min, Max, FX, dFdX, d2FdX2)      \
defineUnivariateEquation2(Name, Type, Min, Max)                                \
addUnivariateEquation2Funcs(Name, Type, FX, dFdX, d2FdX2)

#define createUnivariateEquation3(Name, Type, Min, Max, FX, dFdX, d2FdX2, d3FdX3) \
defineUnivariateEquation3(Name, Type, Min, Max)                                \
addUnivariateEquation3Funcs3(Name, Type, FX, dFdX, d2FdX2, d3FdX3)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
