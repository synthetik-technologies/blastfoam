/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "variableEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::variableEntry::variableEntry()
:
    value_(0.0),
    bounds_(-great, great)
{}



Foam::variableEntry::variableEntry(Istream& is)
:
    value_(0.0),
    bounds_(-great, great)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::variableEntry::~variableEntry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::variableEntry::read(Istream& is)
{
    optimizationEntry::read(is);

    // Read bounds
    is >> bounds_;

    token t(is);
    if (t.isNumber())
    {
        value_ = t.number();
    }
    else
    {
        is.putBack(t);
        value_ = t_->number();
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, variableEntry& entry)
{
    entry.read(is);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const variableEntry& entry)
{
    os  << static_cast<const optimizationEntry&>(entry) << incrIndent
        << indent << "bounds: " << entry.bounds_ << nl
        << indent << "value: " << entry.value_ << decrIndent << endl;
    return os;
}



// ************************************************************************* //
