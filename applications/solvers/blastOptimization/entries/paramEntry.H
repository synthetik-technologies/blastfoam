/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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
    Foam::paramEntry

Description
    optimization entry

SourceFiles
    paramEntry.C

\*---------------------------------------------------------------------------*/
#ifndef paramEntry_H
#define paramEntry_H

#include "optEntry.H"
#include "List.H"

namespace Foam
{

class paramEntry;
Istream& operator>>(Istream& is, paramEntry& entry);
Ostream& operator<<(Ostream& os, const paramEntry& entry);

class paramEntry
:
    public optEntry
{
//- Private data

    //- Value
    token value_;


public:
    //Constructors

        //- Null constructor
        paramEntry();

        //- Construct from Istream
        paramEntry(Istream& is);

    //Destructor
    ~paramEntry();


    // Member functions

        //- Return values
        inline const token& value() const
        {
            return value_;
        }

        //- Access current value
        inline void set()
        {
            (*t_) = value_;
        }

        //- Read
        virtual void read(Istream& is);


    // Operators

        //- Read
        friend Istream& operator>>(Istream& is, paramEntry& entry);

        //- Write
        friend Ostream& operator<<(Ostream& os, const paramEntry& entry);
};


} // End namespace Foam

#endif


// ************************************************************************* //
