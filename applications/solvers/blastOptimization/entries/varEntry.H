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
    Foam::varEntry

Description
    Variable entry

SourceFiles
    varEntry.C

\*---------------------------------------------------------------------------*/
#ifndef varEntry_H
#define varEntry_H

#include "optEntry.H"
#include "Pair.H"

namespace Foam
{

class varEntry;
Istream& operator>>(Istream& is, varEntry& entry);
Ostream& operator<<(Ostream& os, const varEntry& entry);

class varEntry
:
    public optEntry
{
// Private data

    //- Current value
    scalar value_;

    //- Bounds
    Pair<scalar> bounds_;

public:
    //Constructors

        //- Null constructor
        varEntry();

        //- Construct from Istream
        varEntry(Istream& is);

    //Destructor
    ~varEntry();


    // Member functions

        //- Return bounds
        inline const Pair<scalar>& bounds() const
        {
            return bounds_;
        }

        //- Return current value
        inline scalar value() const
        {
            return value_;
        }

        //- Access current value
        inline const token& t()
        {
            return *t_;
        }

        //- Access current value
        inline void set(const scalar& s)
        {
            (*t_) = s;
            value_ = s;
        }

        //- Read
        virtual void read(Istream& is);


    // Operators

        //- Read
        friend Istream& operator>>(Istream& is, varEntry& entry);

        //- Write
        friend Ostream& operator<<(Ostream& os, const varEntry& entry);
};


} // End namespace Foam

#endif


// ************************************************************************* //
