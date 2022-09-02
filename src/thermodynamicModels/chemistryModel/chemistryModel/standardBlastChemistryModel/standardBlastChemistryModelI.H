/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
inline Foam::label
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + density
    return nSpecie_ + 2;
}


template<class BasicThermo, class ThermoType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::RR()
{
    return RR_;
}


template<class BasicThermo, class ThermoType>
inline const Foam::mixtureBlastThermo<BasicThermo, ThermoType>&
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::mixture() const
{
    return mixture_;
}


template<class BasicThermo, class ThermoType>
inline const Foam::PtrList<Foam::Reaction<ThermoType>>&
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::reactions() const
{
    return reactions_;
}


template<class BasicThermo, class ThermoType>
inline const Foam::PtrList<ThermoType>&
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::specieThermos() const
{
    return specieThermos_;
}


template<class BasicThermo, class ThermoType>
inline Foam::label
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::nSpecie() const
{
    return nSpecie_;
}


template<class BasicThermo, class ThermoType>
inline Foam::label
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::nReaction() const
{
    return nReaction_;
}


template<class BasicThermo, class ThermoType>
inline Foam::scalar
 Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::Treact() const
{
    return Treact_;
}


template<class BasicThermo, class ThermoType>
inline Foam::scalar&
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::Treact()
{
    return Treact_;
}


template<class BasicThermo, class ThermoType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::RR(const label i) const
{
    return RR_[i];
}

template<class BasicThermo, class ThermoType>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::RR(const label i)
{
    return RR_[i];
}


// ************************************************************************* //
