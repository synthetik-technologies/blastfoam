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
public:

    virtual word type() const = 0;

    virtual ~indexer()
    {}

    static autoPtr<indexer> New(const List<scalar>&);
    static autoPtr<indexer> New(const word&);
    virtual autoPtr<indexer> clone() const = 0;

    virtual label findIndex(const scalar x, const List<scalar>& xs) const = 0;
};


namespace indexers
{

/*---------------------------------------------------------------------------*\
                        Uniform indexer Declaration
\*---------------------------------------------------------------------------*/\

class uniform
:
    public indexer
{
public:

    virtual ~uniform()
    {}

    virtual word type() const
    {
        return "uniform";
    }

    virtual autoPtr<indexer> clone() const
    {
        return autoPtr<indexer>(new uniform());
    }

    virtual label findIndex(const scalar x, const List<scalar>& xs) const;
};


/*---------------------------------------------------------------------------*\
                        Nonuniform indexer Declaration
\*---------------------------------------------------------------------------*/


class nonuniform
:
    public indexer
{
public:

    virtual ~nonuniform()
    {}

    virtual word type() const
    {
        return "nonuniform";
    }

    virtual autoPtr<indexer> clone() const
    {
        return autoPtr<indexer>(new nonuniform());
    }

    virtual label findIndex(const scalar x, const List<scalar>& xs) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace indexers


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
