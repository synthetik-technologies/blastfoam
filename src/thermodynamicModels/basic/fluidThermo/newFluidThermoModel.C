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

------------------------------------------------------------------------*/

#include "fluidThermoModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermoModel> Foam::fluidThermoModel::NewBasic
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
{
    basicConstructorTable::iterator cstrIter =
        lookupThermo<basicConstructorTable>
        (
            dict,
            basicConstructorTablePtr_
        );

    return cstrIter()(phaseName, p, rho, e, T, dict, master);
}


Foam::autoPtr<Foam::fluidThermoModel> Foam::fluidThermoModel::NewDetonating
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
{
    detonatingConstructorTable::iterator cstrIter =
        lookupThermo<detonatingConstructorTable>
        (
            dict,
            detonatingConstructorTablePtr_
        );

    return cstrIter()(phaseName, p, rho, e, T, dict, master);
}

Foam::autoPtr<Foam::fluidThermoModel> Foam::fluidThermoModel::New
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
{
    word type = dict.lookup("type");
    if (type == "basic")
    {
        return NewBasic(phaseName, p, rho, e, T, dict, master);
    }
    else if (type == "detonating")
    {
        return NewDetonating(phaseName, p, rho, e, T, dict, master);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown fluidThermo type " << type << nl << nl
            << "Valid " << fluidThermoModel::typeName << " types are:" << nl
            << "basic" << nl
            << "detonating" << nl
            << abort(FatalError);
    }
    return autoPtr<fluidThermoModel>();
}

// ************************************************************************* //
