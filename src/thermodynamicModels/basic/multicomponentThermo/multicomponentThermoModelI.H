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

template<class BasicThermo>
inline const Foam::speciesTable&
Foam::multicomponentThermoModel<BasicThermo>::species() const
{
    return species_;
}


template<class BasicThermo>
inline bool Foam::multicomponentThermoModel<BasicThermo>::contains
(
    const word& specieName
) const
{
    return species_.contains(specieName);
}


template<class BasicThermo>
inline bool
Foam::multicomponentThermoModel<BasicThermo>::active(label speciei) const
{
    return active_[speciei];
}


template<class BasicThermo>
inline const Foam::List<bool>&
Foam::multicomponentThermoModel<BasicThermo>::active() const
{
    return active_;
}


template<class BasicThermo>
inline void Foam::multicomponentThermoModel<BasicThermo>::setActive(label speciei)
{
    active_[speciei] = true;
}


template<class BasicThermo>
inline void Foam::multicomponentThermoModel<BasicThermo>::setInactive(label speciei)
{
    active_[speciei] = false;
}


template<class BasicThermo>
inline Foam::PtrList<Foam::volScalarField>&
Foam::multicomponentThermoModel<BasicThermo>::Y()
{
    return Ys_;
}


template<class BasicThermo>
inline const Foam::PtrList<Foam::volScalarField>&
Foam::multicomponentThermoModel<BasicThermo>::Y() const
{
    return Ys_;
}


template<class BasicThermo>
inline Foam::volScalarField&
Foam::multicomponentThermoModel<BasicThermo>::Y(const label i)
{
    return Ys_[i];
}


template<class BasicThermo>
inline const Foam::volScalarField&
Foam::multicomponentThermoModel<BasicThermo>::Y
(
    const label i
) const
{
    return Ys_[i];
}


template<class BasicThermo>
inline Foam::volScalarField&
Foam::multicomponentThermoModel<BasicThermo>::Y
(
    const word& specieName
)
{
    return Ys_[species_[specieName]];
}


template<class BasicThermo>
inline const Foam::volScalarField&
Foam::multicomponentThermoModel<BasicThermo>::Y
(
    const word& specieName
) const
{
    return Ys_[species_[specieName]];
}


// ************************************************************************* //

