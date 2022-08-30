/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "objectMotionState.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::objectMotionState::write(dictionary& dict) const
{
    dict.add("centreOfRotation", centreOfRotation_);
    dict.add("orientation", Q_); // Handled outside
    dict.add("velocity", v_);
    dict.add("acceleration", a_);
    dict.add("angularMomentum", pi_);
    dict.add("torque", tau_);
}


void Foam::objectMotionState::write(Ostream& os) const
{
    writeEntry(os, "centreOfRotation", centreOfRotation_);
    writeEntry(os, "orientation", Q_); // Handled outside
    writeEntry(os, "velocity", v_);
    writeEntry(os, "acceleration", a_);
    writeEntry(os, "angularMomentum", pi_);
    writeEntry(os, "torque", tau_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    objectMotionState& state
)
{
    is  >> state.centreOfRotation_
        >> state.Q_ // Handled outside
        >> state.v_
        >> state.a_
        >> state.pi_
        >> state.tau_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::objectMotionState&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const objectMotionState& state
)
{
    os  << token::SPACE << state.centreOfRotation()
        << token::SPACE << state.Q() // Handled outside
        << token::SPACE << state.v()
        << token::SPACE << state.a()
        << token::SPACE << state.pi()
        << token::SPACE << state.tau();

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::objectMotionState&)"
    );

    return os;
}


// ************************************************************************* //
