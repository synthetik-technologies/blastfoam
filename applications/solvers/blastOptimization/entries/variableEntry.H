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
    Foam::variableEntry

Description
    Variable entry

SourceFiles
    variableEntry.C

\*---------------------------------------------------------------------------*/
#ifndef variableEntry_H
#define variableEntry_H

#include "optimizationEntry.H"
#include "Pair.H"

namespace Foam
{

class variableEntry;
Istream& operator>>(Istream& is, variableEntry& entry);
Ostream& operator<<(Ostream& os, const variableEntry& entry);

class variableEntry
:
    public optimizationEntry
{
// Private data

    //- Current value
    scalar value_;

    //- Bounds
    Pair<scalar> bounds_;

public:
    //Constructors

        //- Null constructor
        variableEntry();

        //- Construct from Istream
        variableEntry(Istream& is);

    //Destructor
    ~variableEntry();


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
        friend Istream& operator>>(Istream& is, variableEntry& entry);

        //- Write
        friend Ostream& operator<<(Ostream& os, const variableEntry& entry);
};


} // End namespace Foam

#endif


// ************************************************************************* //
