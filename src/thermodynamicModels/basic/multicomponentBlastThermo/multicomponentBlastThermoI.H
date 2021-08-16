/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

\*---------------------------------------------------------------------------*/

inline bool
Foam::multicomponentBlastThermo::active(label speciei) const
{
    return active_[speciei];
}


inline const Foam::List<bool>&
Foam::multicomponentBlastThermo::active() const
{
    return active_;
}


inline void Foam::multicomponentBlastThermo::setActive(label speciei)
{
    active_[speciei] = true;
}


inline void Foam::multicomponentBlastThermo::setInactive(label speciei)
{
    active_[speciei] = false;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::multicomponentBlastThermo::Ys()
{
    return Ys_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::multicomponentBlastThermo::Ys() const
{
    return Ys_;
}


inline Foam::volScalarField&
Foam::multicomponentBlastThermo::Y(const label i)
{
    return Ys_[i];
}


inline const Foam::volScalarField& Foam::multicomponentBlastThermo::Y
(
    const label i
) const
{
    return Ys_[i];
}


inline Foam::volScalarField& Foam::multicomponentBlastThermo::Y
(
    const word& specieName
)
{
    return Ys_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::multicomponentBlastThermo::Y
(
    const word& specieName
) const
{
    return Ys_[species_[specieName]];
}


// ************************************************************************* //

