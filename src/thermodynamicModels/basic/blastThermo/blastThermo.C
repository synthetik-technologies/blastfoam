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

#include "blastThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(blastThermo, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastThermo::blastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    timeIntegrationSystem
    (
        IOobject::groupName("blastThermo", phaseName),
        mesh
    ),
    basicThermo::implementation(mesh, dict, phaseName),
    phaseName_(phaseName),
    e_
    (
        IOobject
        (
            basicThermo::phasePropertyName
            (
                "e", phaseName
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass, 0.0),
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    ),
    Cp_
    (
        IOobject
        (
            IOobject::groupName("Cp", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            IOobject::groupName("Cv", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    TLow_
    (
        dict.lookupOrDefault("limitT", true)
      ? dict.lookupOrDefault<scalar>("TLow", 0.0)
      : -great
    ),
    residualAlpha_("residualAlpha", dimless, 0.0),
    residualRho_("residualRho", dimDensity, 0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blastThermo::~blastThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::blastThermo::read()
{
    this->residualRho_.read(*this);
    this->residualAlpha_.read(*this);
    return true;
}


Foam::UIndirectList<Foam::scalar> Foam::blastThermo::cellSetScalarList
(
    const volScalarField& psi,
    const labelList& cells
)
{
    return UIndirectList<scalar>(psi, cells);
}


Foam::word Foam::blastThermo::readThermoType(const dictionary& dict)
{
    return word
    (
        word(dict.lookup("transport")) + '<'
      + word(dict.lookup("thermo")) + '<'
      + word(dict.lookup("equationOfState"))  + '<'
      + word("specie") + ">>>"
    );
}

Foam::wordList Foam::blastThermo::splitThermoName
(
    const word& thermoName,
    const int nCmpt
)
{
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = std::min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");

            // If the number of number of components in the name
            // is greater than nCmpt return an empty list
            if (i == nCmpt)
            {
                return wordList();
            }
        }
        beg = end + 1;
    }

    // If the number of number of components in the name is not equal to nCmpt
    // return an empty list
    if (i + 1 != nCmpt)
    {
        return wordList();
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i].replaceAll(">","");
    }

    return cmpts;
}


Foam::wordList Foam::blastThermo::splitThermoName
(
    const word& thermoName
)
{
    wordList cmpts;
    string::size_type beg=0, end=0, endb=0, endc=0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
          end = std::min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            word newStr = thermoName.substr(beg, end-beg);
            newStr.replaceAll(">","");
            cmpts.append(newStr);
        }
        beg = end + 1;
    }
    if (beg < thermoName.size())
    {
        word newStr = thermoName.substr(beg, string::npos);
        newStr.replaceAll(">","");
        cmpts.append(newStr);
    }

    wordList cmptsFinal(6);
    if (cmpts[0] == "detonating")
    {
        cmptsFinal[0] = cmpts[0];
        cmptsFinal[1] = cmpts[1];
        cmptsFinal[2] = cmpts[2] + '/' + cmpts[6];
        cmptsFinal[3] = cmpts[3] + '/' + cmpts[7];
        cmptsFinal[4] = cmpts[4] + '/' + cmpts[8];
        cmptsFinal[5] = cmpts[5] + '/' + cmpts[9];
    }
    else
    {
        cmptsFinal[0] = cmpts[0];
        cmptsFinal[1] = cmpts[1];
        cmptsFinal[2] = cmpts[2];
        cmptsFinal[3] = cmpts[3];
        cmptsFinal[4] = cmpts[4];
        cmptsFinal[5] = cmpts[5];
    }

    return cmptsFinal;
}


Foam::volScalarField& Foam::blastThermo::lookupOrConstruct
(
    const fvMesh& mesh,
    const word& name,
    const IOobject::readOption rOpt,
    const IOobject::writeOption wOpt,
    const dimensionSet& dims,
    const bool allowNoGroup
)
{
    const word baseName(IOobject::member(name));
    if (!mesh.foundObject<volScalarField>(name))
    {
        volScalarField* fPtr = nullptr;
        IOobject io
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            wOpt
        );
        IOobject baseIo
        (
            baseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            wOpt
        );

        if (io.typeHeaderOk<volScalarField>(true))
        {
            io.readOpt() = IOobject::MUST_READ;
            fPtr =
                new volScalarField
                (
                    io,
                    mesh
                );
        }
        else if (mesh.foundObject<volScalarField>(baseName) && allowNoGroup)
        {
            const volScalarField& baseField =
                mesh.lookupObjectRef<volScalarField>(baseName);
            fPtr =
                new volScalarField
                (
                    io,
                    baseField,
                    baseField.boundaryField()
                );
        }
        else if (baseIo.typeHeaderOk<volScalarField>(true) && allowNoGroup)
        {
            baseIo.readOpt() = IOobject::MUST_READ;
            fPtr =
                new volScalarField
                (
                    baseIo,
                    mesh
                );

            // Rename to the desired name
            fPtr->rename(name);
        }
        else if (rOpt != IOobject::MUST_READ)
        {
            fPtr =
                new volScalarField
                (
                    io,
                    mesh,
                    dimensionedScalar("0", dims, Zero)
                );
        }
        else if (allowNoGroup)
        {
            FatalErrorInFunction
                << name << " or " << baseName << " is required" << endl
                << abort(FatalError);
        }
        else
        {
            FatalErrorInFunction
                << name << " is required" << endl
                << abort(FatalError);
        }

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.lookupObjectRef<volScalarField>(name);
}

const Foam::volScalarField& Foam::blastThermo::he() const
{
    return e_;
}


Foam::volScalarField& Foam::blastThermo::he()
{
    return e_;
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::scalar Foam::blastThermo::cellrho(const label celli) const
{
    return rho_[celli];
}


Foam::volScalarField& Foam::blastThermo::rho()
{
    return rho_;
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::rho0() const
{
    return rho_.oldTime();
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::gamma() const
{
    return volScalarField::New("gamma", Cp_/Cv_);
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::gamma
(
    const scalarField& T,
    const label patchi
) const
{
    return this->Cp(T, patchi)/this->Cv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::kappa() const
{
    return volScalarField::New("kappa", Cp_*this->alpha_);
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::kappa
(
    const label patchi
) const
{
    return
        this->Cp(this->T_.boundaryField()[patchi], patchi)
       *this->alpha_.boundaryField()[patchi];
}


Foam::scalar Foam::blastThermo::cellkappa(const label celli) const
{
    return this->cellCp(this->T_[celli], celli)*this->alpha_[celli];
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::alphahe() const
{
    return volScalarField::New
    (
        "alphahe",
        this->gamma()*this->alpha_
    );
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::alphahe
(
    const label patchi
) const
{
    return
        this->gamma(this->T_.boundaryField()[patchi], patchi)
       *this->alpha_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return Cp_*(this->alpha_ + alphat);
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->Cp(this->T_.boundaryField()[patchi], patchi)
       *(this->alpha_.boundaryField()[patchi] + alphat);
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return this->gamma()*(this->alpha_ + alphat);
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->gamma(this->T_.boundaryField()[patchi], patchi)
       *(this->alpha_.boundaryField()[patchi] + alphat);
}

// ************************************************************************* //
