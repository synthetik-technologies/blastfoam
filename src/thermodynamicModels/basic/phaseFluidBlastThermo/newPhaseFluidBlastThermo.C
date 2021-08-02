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

#include "phaseFluidBlastThermo.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseFluidBlastThermo> Foam::phaseFluidBlastThermo::NewBasic
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
{
    basicConstructorTable::iterator cstrIter =
        lookupThermo<basicConstructorTable>
        (
            dict,
            basicConstructorTablePtr_
        );

    return cstrIter()(phaseName, p, rho, e, T, dict, master, masterName);
}


Foam::autoPtr<Foam::phaseFluidBlastThermo> Foam::phaseFluidBlastThermo::NewDetonating
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
{
    detonatingConstructorTable::iterator cstrIter =
        lookupThermo<detonatingConstructorTable>
        (
            dict,
            detonatingConstructorTablePtr_
        );

    return cstrIter()(phaseName, p, rho, e, T, dict, master, masterName);
}


Foam::autoPtr<Foam::phaseFluidBlastThermo> Foam::phaseFluidBlastThermo::NewMulticomponent
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
{
    multicomponentConstructorTable::iterator cstrIter =
        lookupThermo<multicomponentConstructorTable>
        (
            dict,
            multicomponentConstructorTablePtr_
        );

    return cstrIter()(phaseName, p, rho, e, T, dict, master, masterName);
}


Foam::autoPtr<Foam::phaseFluidBlastThermo> Foam::phaseFluidBlastThermo::NewReacting
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
{
    reactingConstructorTable::iterator cstrIter =
        lookupThermo<multicomponentConstructorTable>
        (
            dict,
            reactingConstructorTablePtr_
        );

    return cstrIter()(phaseName, p, rho, e, T, dict, master, masterName);
}


Foam::autoPtr<Foam::phaseFluidBlastThermo> Foam::phaseFluidBlastThermo::New
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
{
    const word type(dict.lookup<word>("type"));
    if (type == "basic")
    {
        return NewBasic(phaseName, p, rho, e, T, dict, master, masterName);
    }
    else if (type == "detonating")
    {
        return NewDetonating(phaseName, p, rho, e, T, dict, master, masterName);
    }
    else if (type == "multicomponent")
    {
        return NewMulticomponent
        (
            phaseName, p, rho, e, T, dict, master, masterName
        );
    }
//     else if (type == "reacting")
//     {
//         return NewReacting
//         (
//             phaseName, p, rho, e, T, dict, master, masterName
//         );
//     }
    else
    {
        FatalErrorInFunction
            << "Unknown fluidThermo type " << type << nl << nl
            << "Valid " << phaseFluidBlastThermo::typeName << " types are:" << nl
            << "basic" << nl
            << "detonating" << nl
            << "multicomponent" << nl
//             << "reacting" << nl
            << abort(FatalError);
    }
    return autoPtr<phaseFluidBlastThermo>();
}

// ************************************************************************* //
