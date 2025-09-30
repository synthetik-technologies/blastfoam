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

    //- Reference to indexed data
    const List<scalar>& xs_;


public:

    //- Type information
    TypeName("indexer");

    //- Constructor
    indexer(const List<scalar>& xs);

    //- Clone the indexer
    virtual autoPtr<indexer> clone() const = 0;


    //- Destructor
    virtual ~indexer()
    {}

    // Static Member Functions

        //- Selector based on the list
        static autoPtr<indexer> New
        (
            const List<scalar>&
        );

        //- Select indexer based on the type
        static autoPtr<indexer> New
        (
            const word&,
            const List<scalar>&
        );


    // Public Member Functions

        //- Find  the lower index given a position
        virtual label findIndex(const scalar x) const = 0;

        //- Validate the indexer
        virtual bool validate(bool fail = true) const;
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

    //- Type information
    TypeName("null");

    //- Constructor
    null(const List<scalar>& xs)
    :
        indexer(xs)
    {}

    //- Clone the indexer
    virtual autoPtr<indexer> clone() const
    {
        return autoPtr<indexer>(new null(xs_));
    }


    //- Destructor
    virtual ~null()
    {}


    // Public Member Functions

        //- Find  the lower index given a position
        inline label findIndex(const scalar x) const
        {
            return 0;
        }

        //- Validate the indexer
        virtual bool validate(bool fail = true) const
        {
            return true;
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

    //- Type information
    TypeName("uniform");

    //- Constructor
    uniform(const List<scalar>& xs)
    :
        indexer(xs)
    {}

    //- Clone the indexer
    virtual autoPtr<indexer> clone() const
    {
        return autoPtr<indexer>(new uniform(xs_));
    }


    //- Destructor
    virtual ~uniform()
    {}


    // Public Member Functions

        //- Find  the lower index given a position
        virtual label findIndex(const scalar x) const;
};


/*---------------------------------------------------------------------------*\
                        Nonuniform indexer Declaration
\*---------------------------------------------------------------------------*/


class nonuniform
:
    public indexer
{
// Private member data

    //- Cached index
    mutable label i_;


public:

    //- Type information
    TypeName("nonuniform");

    //- Constructor
    nonuniform(const List<scalar>& xs)
    :
        indexer(xs)
    {}

    //- Clone the indexer
    virtual autoPtr<indexer> clone() const
    {
        return autoPtr<indexer>(new nonuniform(xs_));
    }


    //- Destructor
    virtual ~nonuniform()
    {}


    // Public Member Functions

        //- Find  the lower index given a position
        virtual label findIndex(const scalar x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace indexers


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
