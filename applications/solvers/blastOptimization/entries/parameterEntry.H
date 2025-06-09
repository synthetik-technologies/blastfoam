/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023-2025
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
    Foam::parameterEntry

Description
    optimization entry

SourceFiles
    parameterEntry.C

\*---------------------------------------------------------------------------*/
#ifndef parameterEntry_H
#define parameterEntry_H

#include "optimizationEntry.H"
#include "List.H"

namespace Foam
{

class parameterEntry;
Istream& operator>>(Istream& is, parameterEntry& entry);
Ostream& operator<<(Ostream& os, const parameterEntry& entry);

class parameterEntry
:
    public optimizationEntry
{
//- Private data

    //- Value
    token value_;


public:
    //Constructors

        //- Null constructor
        parameterEntry();

        //- Construct from Istream
        parameterEntry(Istream& is);

    //Destructor
    ~parameterEntry();


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
        friend Istream& operator>>(Istream& is, parameterEntry& entry);

        //- Write
        friend Ostream& operator<<(Ostream& os, const parameterEntry& entry);
};


} // End namespace Foam

#endif


// ************************************************************************* //
