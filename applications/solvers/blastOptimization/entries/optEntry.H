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
    Foam::optEntry

Description
    optimization entry

SourceFiles
    optEntry.C

\*---------------------------------------------------------------------------*/
#ifndef optEntry_H
#define optEntry_H

#include "HashTable.H"
#include "dictionary.H"

namespace Foam
{

class optEntry;
Istream& operator>>(Istream& is, optEntry& entry);
Ostream& operator<<(Ostream& os, const optEntry& entry);

class optEntry
{
public:

    static HashTable<dictionary> dictionaries;

protected:

    //- Name of file this entry is read from
    fileName file_;

    //- Path to the token
    fileName path_;

    //- Index (if in a list)
    label index_;

    //- Pointer to the actual token for easy access
    token* t_;


    // Private member functions

        //- Find and return the token corresponding to the entry
        static token* findTokenEntry
        (
            const dictionary& dict,
            const fileName& path,
            label index
        );

        //- Return the dictionary corresponding to the given path
        static fileName findDictFromPath(fileName& path);


public:
    //Constructors

        //- Null constructor
        optEntry();

        //- Construct from Istream
        optEntry(Istream& is);

    //Destructor
    virtual ~optEntry();


    // Member functions

        //- Return file
        inline const fileName& file() const
        {
            return file_;
        }

        //- Return path
        inline const fileName& path() const
        {
            return path_;
        }

        //- Access current value
        inline const token& t()
        {
            return *t_;
        }

        //- Read
        virtual void read(Istream& is);


    // Operators

        virtual bool operator!=(const optEntry& e) const
        {
            return t_ != e.t_;
        }
        virtual bool operator==(const optEntry& e) const
        {
            return t_ == e.t_;
        }

        //- Read
        friend Istream& operator>>(Istream& is, optEntry& entry);

        //- Write
        friend Ostream& operator<<(Ostream& os, const optEntry& entry);
};


} // End namespace Foam

#endif


// ************************************************************************* //
