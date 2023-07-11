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
#include "runTimeSelectionTables.H"
#include "dictionary.H"
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

    enum RealType
    {
        UNKNOWN,
        REAL,
        MODIFIED
    };

protected:
    // Protected data

        //- Optionally read "isReal" flag"
        RealType isReal_;

public:

    TypeName(pTraits<Type>::typeName + "Modifier");

    Modifier()
    :
        isReal_(UNKNOWN)
    {}
    Modifier(const dictionary& dict)
    :
        isReal_(UNKNOWN)
    {}

    virtual ~Modifier() {}

    declareRunTimeSelectionTable
    (
        autoPtr,
        Modifier,
        null,
        (),
        ()
    );
    declareRunTimeSelectionTable
    (
        autoPtr,
        Modifier,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    static autoPtr<Modifier<Type>> New(const word& mod);
    static autoPtr<Modifier<Type>> New(const word& mod, const dictionary& dict);

    virtual autoPtr<Modifier<Type>> clone() const = 0;

    virtual Type operator()(const Type&) const = 0;
    virtual void Mod(UList<Type>&) const;

    virtual Type inv(const Type&) const = 0;
    virtual void Inv(UList<Type>&) const;

    virtual bool needMod() const
    {
        return true;
    }

    virtual bool isReal() const
    {
        if (isReal_ == UNKNOWN)
        {
            FatalErrorInFunction
                << "Unable to determine if the read data is real" << endl
                << abort(FatalError);
        }
        return isReal_ == REAL;
    }
    virtual void setReal()
    {
        isReal_ = REAL;
    }
    virtual void setMod()
    {
        isReal_ = MODIFIED;
    }
    virtual void readReal(const dictionary& dict, const word& entry)
    {
        isReal_ = dict.template lookup<bool>(entry) ? REAL : MODIFIED;
    }
};

#define defineModifierType(Type)                                               \
typedef Modifier<Type> Type##Modifier;                                         \
defineTemplateTypeNameAndDebug(Type##Modifier, 0);                             \
defineTemplateRunTimeSelectionTable(Type##Modifier, null);                     \
defineTemplateRunTimeSelectionTable(Type##Modifier, dictionary);

#define addDictModifier(Mod, Type)                                             \
defineTypeNameAndDebug(Mod, 0);                                                \
addToRunTimeSelectionTable(Type##Modifier, Mod, dictionary);

#define addDictModifierType(Mod, Type)                                         \
typedef Mod<Type> Type##Mod;                                                   \
defineNamedTemplateTypeNameAndDebug(Type##Mod, 0);                             \
addToRunTimeSelectionTable(Type##Modifier, Type##Mod, dictionary);

#define addModifier(Mod, Type)                                                 \
defineTypeNameAndDebug(Mod, 0);                                                \
addToRunTimeSelectionTable(Type##Modifier, Mod, null);                         \
addToRunTimeSelectionTable(Type##Modifier, Mod, dictionary);

#define addModifierType(Mod, Type)                                             \
typedef Mod<Type> Type##Mod;                                                   \
defineNamedTemplateTypeNameAndDebug(Type##Mod, 0);                             \
addToRunTimeSelectionTable(Type##Modifier, Type##Mod, null);                   \
addToRunTimeSelectionTable(Type##Modifier, Type##Mod, dictionary);

#define addDictModifiers(Mod)                                                  \
addDictModifierType(Mod, scalar)                                               \
addDictModifierType(Mod, vector)                                               \
addDictModifierType(Mod, symmTensor)                                           \
addDictModifierType(Mod, sphericalTensor)                                      \
addDictModifierType(Mod, tensor)

#define addModifiers(Mod)                                                      \
addModifierType(Mod, scalar)                                                   \
addModifierType(Mod, vector)                                                   \
addModifierType(Mod, symmTensor)                                               \
addModifierType(Mod, sphericalTensor)                                          \
addModifierType(Mod, tensor)

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

    TypeName("none");

    None()
    {
        this->isReal_ = Modifier<Type>::REAL;
    }
    None(const dictionary& dict)
    {
        this->isReal_ = Modifier<Type>::REAL;
    }



    virtual autoPtr<Modifier<Type>> clone() const
    {
        return autoPtr<Modifier<Type>>(new None<Type>());
    }

    inline Type operator()(const Type& x) const
    {
        return x;
    }
    virtual void Mod(UList<Type>&) const
    {}

    inline Type inv(const Type& x) const
    {
        return x;
    }
    virtual void Inv(UList<Type>&) const
    {}

    inline bool needMod() const
    {
        return false;
    }

    inline bool isReal() const
    {
        return true;
    }
    inline void setReal()
    {}
    inline void setMod()
    {}
    virtual void readReal(const dictionary& dict, const word& entry)
    {}
};


/*---------------------------------------------------------------------------*\
                        Scalar modifiers Declaration
\*---------------------------------------------------------------------------*/


class Ln
:
    public Modifier<scalar>
{
public:

    TypeName("ln");

    using Modifier<scalar>::Modifier;

    virtual autoPtr<Modifier<scalar>> clone() const
    {
        return autoPtr<Modifier<scalar>>(new Ln());
    }

    inline scalar operator()(const scalar& x) const
    {
        return log(max(x, small));
    }

    inline scalar inv(const scalar& x) const
    {
        return exp(x);
    }
};


class Exp
:
    public Modifier<scalar>
{
public:

    TypeName("exp");

    using Modifier<scalar>::Modifier;

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
        return log(max(x, small));
    }
};


class Log10
:
    public Modifier<scalar>
{
public:

    TypeName("log10");

    using Modifier<scalar>::Modifier;

    virtual autoPtr<Modifier<scalar>> clone() const
    {
        return autoPtr<Modifier<scalar>>(new Log10());
    }

    inline scalar operator()(const scalar& x) const
    {
        return log10(max(x, small));
    }

    inline scalar inv(const scalar& x) const
    {
        return pow(10.0, x);
    }
};


class Pow10
:
    public Modifier<scalar>
{
public:

    TypeName("pow10");

    using Modifier<scalar>::Modifier;

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
        return log10(max(x, small));
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace modifiers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ModifierTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

