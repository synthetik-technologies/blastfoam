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
    Foam::indexer
    Foam::indexers::uniform
    Foam::indexers::nonuniform

Description
    Collection of indexers

SourceFiles
    indexing.C

\*---------------------------------------------------------------------------*/

#ifndef indexing_H
#define indexing_H

#include "scalar.H"
#include "List.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class indexer Declaration
\*---------------------------------------------------------------------------*/


class indexer
{
protected:
    const List<scalar>& xs_;

public:

    virtual word type() const = 0;

    indexer(const List<scalar>& xs)
    :
        xs_(xs)
    {}

    virtual ~indexer()
    {}

    static autoPtr<indexer> New(const List<scalar>&);
    static autoPtr<indexer> New(const word&, const List<scalar>&);
    virtual autoPtr<indexer> clone() const = 0;

    virtual label findIndex(const scalar x) const = 0;
};


namespace indexers
{

/*---------------------------------------------------------------------------*\
                        null indexer Declaration
\*---------------------------------------------------------------------------*/\

class null
:
    public indexer
{
public:

    null(const List<scalar>& xs)
    :
        indexer(xs)
    {}
    virtual ~null()
    {}

    virtual word type() const
    {
        return "null";
    }

    virtual autoPtr<indexer> clone() const
    {
        return autoPtr<indexer>(new null(xs_));
    }

    virtual label findIndex(const scalar x) const
    {
        return 0;
    }
};


/*---------------------------------------------------------------------------*\
                        Uniform indexer Declaration
\*---------------------------------------------------------------------------*/\

class uniform
:
    public indexer
{
public:

    uniform(const List<scalar>& xs)
    :
        indexer(xs)
    {}

    virtual ~uniform()
    {}

    virtual word type() const
    {
        return "uniform";
    }

    virtual autoPtr<indexer> clone() const
    {
        return autoPtr<indexer>(new uniform(xs_));
    }

    virtual label findIndex(const scalar x) const;
};


/*---------------------------------------------------------------------------*\
                        Nonuniform indexer Declaration
\*---------------------------------------------------------------------------*/


class nonuniform
:
    public indexer
{
public:

    nonuniform(const List<scalar>& xs)
    :
        indexer(xs)
    {}

    virtual ~nonuniform()
    {}

    virtual word type() const
    {
        return "nonuniform";
    }

    virtual autoPtr<indexer> clone() const
    {
        return autoPtr<indexer>(new nonuniform(xs_));
    }

    virtual label findIndex(const scalar xff) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace indexers


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
