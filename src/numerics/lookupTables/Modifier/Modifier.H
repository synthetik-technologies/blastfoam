/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

Class
    Foam::Modifier

Description
    Collection of modifiers

SourceFiles
    Modifiers.C

\*---------------------------------------------------------------------------*/

#ifndef Modifiers_H
#define Modifiers_H

#include "scalar.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Modifier Declaration
\*---------------------------------------------------------------------------*/


template<class Type>
class Modifier
{
public:

    virtual word type() const = 0;

    virtual ~Modifier()
    {}

    static autoPtr<Modifier<Type>> New(const word& mod);
    virtual autoPtr<Modifier<Type>> clone() const = 0;

    virtual Type operator()(const Type&) const = 0;
    virtual Type inv(const Type&) const = 0;

    virtual bool needMod() const = 0;
};

// Specilization of scalar modifier
template<>
autoPtr<Modifier<scalar>> Modifier<scalar>::New(const word& mod);

/*---------------------------------------------------------------------------*\
                        Templated modifiers Declaration
\*---------------------------------------------------------------------------*/

namespace modifiers
{

template<class Type>
class None
:
    public Modifier<Type>
{
public:

    virtual word type() const
    {
        return "none";
    }

    virtual autoPtr<Modifier<Type>> clone() const
    {
        return autoPtr<Modifier<Type>>(new None<Type>());
    }

    inline Type operator()(const Type& x) const
    {
        return x;
    }

    inline Type inv(const Type& x) const
    {
        return x;
    }

    inline bool needMod() const
    {
        return false;
    }
};


/*---------------------------------------------------------------------------*\
                        Scalar modifiers Declaration
\*---------------------------------------------------------------------------*/


class Ln
:
    public Modifier<scalar>
{
public:

    virtual word type() const
    {
        return "ln";
    }

    virtual autoPtr<Modifier<scalar>> clone() const
    {
        return autoPtr<Modifier<scalar>>(new Ln());
    }

    inline scalar operator()(const scalar& x) const
    {
        return log(x);
    }

    inline scalar inv(const scalar& x) const
    {
        return exp(x);
    }

    inline bool needMod() const
    {
        return true;
    }
};


class Exp
:
    public Modifier<scalar>
{
public:

    virtual word type() const
    {
        return "exp";
    }

    virtual autoPtr<Modifier<scalar>> clone() const
    {
        return autoPtr<Modifier<scalar>>(new Exp());
    }

    inline scalar operator()(const scalar& x) const
    {
        return exp(x);
    }

    inline scalar inv(const scalar& x) const
    {
        return log(x);
    }

    inline bool needMod() const
    {
        return true;
    }
};


class Log10
:
    public Modifier<scalar>
{
public:

    virtual word type() const
    {
        return "log10";
    }

    virtual autoPtr<Modifier<scalar>> clone() const
    {
        return autoPtr<Modifier<scalar>>(new Log10());
    }

    inline scalar operator()(const scalar& x) const
    {
        return log10(x);
    }

    inline scalar inv(const scalar& x) const
    {
        return pow(10.0, x);
    }

    inline bool needMod() const
    {
        return true;
    }
};


class Pow10
:
    public Modifier<scalar>
{
public:

    virtual word type() const
    {
        return "pow10";
    }

    virtual autoPtr<Modifier<scalar>> clone() const
    {
        return autoPtr<Modifier<scalar>>(new Pow10());
    }

    inline scalar operator()(const scalar& x) const
    {
        return pow(10.0, x);
    }

    inline scalar inv(const scalar& x) const
    {
        return log10(x);
    }

    inline bool needMod() const
    {
        return true;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace modifiers


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

