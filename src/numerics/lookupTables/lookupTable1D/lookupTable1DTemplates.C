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

#include "lookupTable1D.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
template<class fType>
void Foam::lookupTable1D<Type>::readComponent
(
    const dictionary& parentDict,
    const word& name,
    word& type,
    Field<fType>& values,
    Field<fType>& modValues,
    modFuncType& modFunc,
    modFuncType& invModFunc,
    const bool canRead
)
{
    Switch isReal = true;
    bool canSetMod = true;
    if (parentDict.found(name + "Coeffs"))
    {
        const dictionary& dict(parentDict.subDict(name + "Coeffs"));
        type = dict.lookupOrDefault<word>("mod", "none");
        setMod
        (
            type,
            modFunc,
            invModFunc
        );

        if (dict.found(name))
        {
            if (canRead)
            {
                values = dict.lookup<Field<fType>>(name);
                modValues.resize(values.size());
                isReal = dict.lookupOrDefault<Switch>("isReal", true);
            }
            else
            {
                canSetMod = false;
            }
        }
        else if (dict.found("file"))
        {
            fileName file(dict.lookup("file"));
            read1DTable
            (
                file,
                dict.lookupOrDefault<string>("delim", ","),
                values
            );
            isReal = dict.lookupOrDefault<Switch>("isReal", true);
        }
        else
        {
            label ny = dict.lookup<label>("n");
            fType dy = dict.lookup<fType>("delta");
            fType miny = dict.lookup<fType>("min");

            values.resize(ny);
            forAll(values, j)
            {
                values[j] = miny + dy*j;
            }
            isReal = dict.lookupOrDefault<Switch>("isReal", true);
        }
    }
    else if (parentDict.found(name))
    {
        if (canRead)
        {
            values = parentDict.lookup<Field<fType>>(name);
            modValues.resize(values.size());
            isReal =
                parentDict.lookupOrDefault<Switch>
                (
                    name + "isReal",
                    true
                );
        }
        else
        {
            canSetMod = false;
        }
        type = parentDict.lookupOrDefault<word>(name + "Mod", "none");
        setMod
        (
            type,
            modFunc,
            invModFunc
        );
    }
    else
    {
        FatalErrorInFunction
            << "Either a list of values or " << name << "Coeffs must" << nl
            << "be provided for " << name << endl
            << abort(FatalError);
    }

    if (canSetMod)
    {
        if (!isReal)
        {
            modValues = values;
            forAll(values, i)
            {
                values[i] = invModFunc(values[i]);
            }
        }
        else
        {
            modValues.resize(values.size());
            forAll(values, i)
            {
                modValues[i] = modFunc(values[i]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<class fType>
fType Foam::lookupTable1D<Type>::interpolate
(
    const fType& fm,
    const fType& fp
) const
{
    return f_ == 0 ? fm : (f_ == 1 ? fp : (1.0 - f_)*fm + f_*fp);
}


template<class Type>
template<template<class> class ListType, class fType>
fType Foam::lookupTable1D<Type>::interpolate
(
    const scalar x,
    const ListType<fType>& fs
) const
{
    update(x);
    return
        f_ == 0 ? fs[index_]
      : (
            f_ == 1
          ? fs[index_ + 1]
          : (1.0 - f_)*fs[index_] + f_*fs[index_ + 1]
        );
}

// ************************************************************************* //
