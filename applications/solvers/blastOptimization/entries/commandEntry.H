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
    Foam::commandEntry

Description
    optimization entry

SourceFiles
    commandEntry.C

\*---------------------------------------------------------------------------*/
#ifndef commandEntry_H
#define commandEntry_H

#include "HashTable.H"
#include "dictionary.H"

namespace Foam
{

class commandEntry;
Istream& operator>>(Istream& is, commandEntry& entry);

class commandEntry
:
    public string
{
public:

    //- Name of file this entry is read from
    word command_;

    //- Options
    List<Tuple2<word, string>> options_;

    //- Number of processors
    label nProcs_;


public:
    //Constructors

        //- Null constructor
        commandEntry();

        //- Construct from Istream
        commandEntry(Istream& is);

    //Destructor
    virtual ~commandEntry();


    // Member functions

        //- Return file
        inline const word& command() const
        {
            return command_;
        }

        //- Return path
        inline const List<Tuple2<word, string>>& options() const
        {
            return options_;
        }


    // Operators

        virtual bool operator!=(const commandEntry& e) const
        {
            return
                command_ != e.command_
             && options_ != e.options_;
        }
        virtual bool operator==(const commandEntry& e) const
        {
            return !(*this != e);
        }

        //- Read
        friend Istream& operator>>(Istream& is, commandEntry& entry);
};


} // End namespace Foam

#endif


// ************************************************************************* //
