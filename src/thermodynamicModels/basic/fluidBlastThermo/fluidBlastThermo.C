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

#include "fluidBlastThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidBlastThermo, 0);
    defineRunTimeSelectionTable(fluidBlastThermo, dictionary);
    defineRunTimeSelectionTable(fluidBlastThermo, phase);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidBlastThermo::fluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName,
    const bool requireRho
)
:
    physicalProperties(mesh, phaseName),
    blastThermo(mesh, dict, phaseName),
    p_
    (
        blastThermo::lookupOrConstruct
        (
            mesh,
            basicThermo::phasePropertyName("p", phaseName),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            dimPressure,
            true // pressure is always allowed to read from mixture field
        )
    ),
    mu_
    (
        IOobject
        (
            basicThermo::phasePropertyName("thermo:mu", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    speedOfSound_
    (
        IOobject
        (
            basicThermo::phasePropertyName("speedOfSound", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimVelocity, 0.0)
    )
{
    if (requireRho && !this->rho_.headerOk())
    {
        FatalErrorInFunction
            << this->rho_.name() << " must be proved for single phase simulations"
            << ", i.e. " << this->rho_.path()/this->rho_.name()
            << endl
            << abort(FatalError);
    }
    this->properties().dictionary::operator=(dict);

    // Make sure the physicalProperties file is not reread
    this->properties().checkOut();
    this->properties().readOpt() = IOobject::MUST_READ;
    this->properties().checkIn();
}


void Foam::fluidBlastThermo::initializeFields()
{
    if (!e_.headerOk())
    {
        const word initType
        (
            this->lookupOrDefault<word>("eInitialization", "pRho")
        );

        //- Calculate internal energy if it was not read
        if (initType == "pRho")
        {
            e_ == this->calce(p_);
        }
        else if (initType == "TRho")
        {
            p_ == this->pRhoT();
            e_ == this->he(p_, T_);
        }
        else
        {
            FatalIOErrorInFunction(*this)
                << "Invalid method of internal energy initialization" << nl
                << "Valid methods are:" << nl
                << "    pRho" << nl
                << "    TRho" << nl
                << endl
                << abort(FatalIOError);
        }
    }
    this->correct();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidBlastThermo> Foam::fluidBlastThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& thermoType,
    const word& phaseName
)
{
    if (thermoType == word::null)
    {
        return blastThermo::New<fluidBlastThermo>
        (
            mesh,
            (dict.isDict(phaseName) && phaseName != word::null)
          ? dict.subDict(phaseName)
          : dict.optionalSubDict("mixture"),
            phaseName,
            phaseName
        );
    }

    phaseConstructorTable::iterator cstrIter =
        phaseConstructorTablePtr_->find(thermoType);

    if (cstrIter == phaseConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluidThermo type " << endl
            << phaseConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh, dict, phaseName);
}


Foam::autoPtr<Foam::fluidBlastThermo> Foam::fluidBlastThermo::New
(
    const fvMesh& mesh,
    const word& thermoType,
    const word& phaseName
)
{
    const IOdictionary dict
    (
        physicalProperties::findModelDict(mesh, phaseName)
    );

    if (thermoType == word::null)
    {
        return blastThermo::New<fluidBlastThermo>
        (
            mesh,
            dict.optionalSubDict("mixture"),
            phaseName,
            phaseName
        );
    }

    phaseConstructorTable::iterator cstrIter =
        phaseConstructorTablePtr_->find(thermoType);

    if (cstrIter == phaseConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluidThermo type " << endl
            << phaseConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh, dict, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidBlastThermo::~fluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidBlastThermo::updateRho()
{
    updateRho(p_);
}


Foam::tmp<Foam::volScalarField> Foam::fluidBlastThermo::pRhoT() const
{
    tmp<volScalarField> tp
    (
        volScalarField::New
        (
            IOobject::groupName("p", this->phaseName()),
            this->rho_.mesh(),
            dimensionedScalar(dimPressure, 0.0)
        )
    );
    volScalarField& p = tp.ref();
    forAll(p, celli)
    {
        p[celli] = this->cellpRhoT(celli);
    }

    volScalarField::Boundary& bp = p.boundaryFieldRef();
    forAll(bp, patchi)
    {
        forAll(bp[patchi], facei)
        {
            bp[patchi][facei] = this->patchFacepRhoT(patchi, facei);
        }
    }
    return tp;
}


Foam::volScalarField& Foam::fluidBlastThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::fluidBlastThermo::p() const
{
    return p_;
}


Foam::tmp<Foam::volScalarField> Foam::fluidBlastThermo::renameRho()
{
    rho_.rename
    (
        basicThermo::phasePropertyName
        (
            Foam::typedName<fluidBlastThermo>("rho")
        )
    );
    return rho_;
}


const Foam::volScalarField& Foam::fluidBlastThermo::speedOfSound() const
{
    return speedOfSound_;
}


Foam::volScalarField& Foam::fluidBlastThermo::speedOfSound()
{
    return speedOfSound_;
}


const Foam::volScalarField& Foam::fluidBlastThermo::mu() const
{
    return mu_;
}


Foam::scalar Foam::fluidBlastThermo::cellnu
(
    const label celli
) const
{
    return mu_[celli]/cellrho(celli);
}


// ************************************************************************* //
