/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "equationBase.H"
#include "OStringStream.H"

// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

Foam::string Foam::equationBase::mergeStrings(const List<string>& eqns) const
{
    OStringStream os;
    os << word(eqns[0]);
    for (label i = 1; i < eqns.size(); i++)
    {
        os << nl << word(eqns[i]);
    }
    return os.str();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationBase::equationBase()
:
    eqnString_(string::null)
{}


Foam::equationBase::equationBase(const string& eqnString)
:
    eqnString_(eqnString)
{}


Foam::equationBase::equationBase(const List<string>& eqnStrings)
:
    eqnString_(mergeStrings(eqnStrings))
{}


Foam::equationBase::equationBase(const dictionary& dict)
:
    eqnString_(dict.lookupOrDefault<string>("eqnString", string::null))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationBase::~equationBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
