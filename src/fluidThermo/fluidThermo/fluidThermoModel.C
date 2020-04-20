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

#include "fluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "blastFixedEnergyFvPatchScalarField.H"
#include "blastGradientEnergyFvPatchScalarField.H"
#include "blastGradientEnergyCalculatedTemperatureFvPatchScalarField.H"
#include "blastMixedEnergyFvPatchScalarField.H"
#include "blastMixedEnergyCalculatedTemperatureFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "blastEnergyJumpFvPatchScalarField.H"
#include "blastEnergyJumpAMIFvPatchScalarField.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermoModel, 0);
    defineRunTimeSelectionTable(fluidThermoModel, basic);
    defineRunTimeSelectionTable(fluidThermoModel, detonating);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermoModel::fluidThermoModel
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName("fluidThermo", phaseName),
            e.mesh().time().timeName(),
            e.mesh()
        )
    ),
    master_(master),
    name_(phaseName),
    p_(p),
    rho_(rho),
    e_(e),
    T_(T),
    mu_
    (
        IOobject
        (
            IOobject::groupName("mu", name_),
            p_.mesh().time().timeName(),
            p_.mesh()
        ),
        p_.mesh(),
        dimensionedScalar(dimDynamicViscosity, 0.0),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    alpha_
    (
        IOobject
        (
            IOobject::groupName("thermo:alpha", name_),
            p_.mesh().time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    residualAlpha_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualAlpha",
            dict,
            dimless,
            0.0
        )
    ),
    residualRho_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualRho",
            dict,
            dimDensity,
            0.0
        )
    ),
    limit_(dict.lookupOrDefault("limit", true))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermoModel::~fluidThermoModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::fluidThermoModel::readThermoType(const dictionary& dict)
{
    return word
    (
        word(dict.lookup("transport")) + '<'
      + word(dict.lookup("thermo")) + '<'
      + word(dict.lookup("equationOfState"))  + '<'
      + word("specie") + ">>>"
    );
}


Foam::wordList Foam::fluidThermoModel::splitThermoName
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

    wordList cmptsFinal(5);
    if (cmpts[0] == "detonating")
    {
        cmptsFinal[0] = "y";
        cmptsFinal[1] = cmpts[1] + '/' + cmpts[5];
        cmptsFinal[2] = cmpts[2] + '/' + cmpts[6];
        cmptsFinal[3] = cmpts[3] + '/' + cmpts[7];
        cmptsFinal[4] = cmpts[4] + '/' + cmpts[8];
    }
    else
    {
        cmptsFinal[0] = "n";
        cmptsFinal[1] = cmpts[0];
        cmptsFinal[2] = cmpts[1];
        cmptsFinal[3] = cmpts[2];
        cmptsFinal[4] = cmpts[3];
    }

    return cmptsFinal;
}

Foam::wordList
Foam::fluidThermoModel::eBoundaryTypes(const volScalarField& T)
{
    wordList ebf = T.boundaryField().types();

    forAll(ebf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(T.boundaryField()[patchi]))
        {
            ebf[patchi] = blastFixedEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<fixedGradientFvPatchScalarField>(T.boundaryField()[patchi])
         || isA<blastGradientEnergyCalculatedTemperatureFvPatchScalarField>
            (
                T.boundaryField()[patchi]
            )
        )
        {
            ebf[patchi] = blastGradientEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<mixedFvPatchScalarField>(T.boundaryField()[patchi])
         || isA<blastMixedEnergyCalculatedTemperatureFvPatchScalarField>
            (
                T.boundaryField()[patchi]
            )
        )
        {
            ebf[patchi] = blastMixedEnergyFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpFvPatchScalarField>(T.boundaryField()[patchi]))
        {
            ebf[patchi] = blastEnergyJumpFvPatchScalarField::typeName;
        }
        else if
        (
            isA<fixedJumpAMIFvPatchScalarField>(T.boundaryField()[patchi])
        )
        {
            ebf[patchi] = blastEnergyJumpAMIFvPatchScalarField::typeName;
        }
    }

    return ebf;
}


Foam::wordList
Foam::fluidThermoModel::eBoundaryBaseTypes(const volScalarField& T)
{
    const volScalarField::Boundary& tbf = T.boundaryField();

    wordList ebt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
        if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpFvPatchScalarField&>(tbf[patchi]);

            ebt[patchi] = pf.interfaceFieldType();
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpAMIFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpAMIFvPatchScalarField&>
                (
                    tbf[patchi]
                );

            ebt[patchi] = pf.interfaceFieldType();
        }
    }

    return ebt;
}



void Foam::fluidThermoModel::eBoundaryCorrection()
{
    volScalarField::Boundary& eBf = e_.boundaryFieldRef();

    forAll(eBf, patchi)
    {
        if (isA<blastGradientEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<blastGradientEnergyFvPatchScalarField>(eBf[patchi]).gradient()
                = eBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<blastMixedEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<blastMixedEnergyFvPatchScalarField>(eBf[patchi]).refGrad()
                = eBf[patchi].fvPatchField::snGrad();
        }
    }
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::nu() const
{
    return mu_/rho_;
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::nu(const label patchi) const
{
    return mu(patchi)/rho_.boundaryField()[patchi];
}


Foam::scalar Foam::fluidThermoModel::nui(const label celli) const
{
    return mu_[celli]/rho_[celli];
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::alpha() const
{
    return alpha_;
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::alphaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New
    (
        "alphaEff",
        CpByCv()*(alpha_ + alphat)
    );
}


Foam::tmp<Foam::scalarField> Foam::fluidThermoModel::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return CpByCv(patchi)*(alpha(patchi) + alphat);
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::alphahe() const
{
    return volScalarField::New
    (
        "alphahe",
        CpByCv()*alpha_
    );
}

Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::alphahe(const label patchi) const
{
    return CpByCv(patchi)*alpha(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::kappa() const
{
    return volScalarField::New
    (
        "kappa",
        Cp()*alpha_
    );
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::kappa(const label patchi) const
{
    return Cp(patchi)*alpha_.boundaryField()[patchi];
}


Foam::scalar Foam::fluidThermoModel::kappai(const label celli) const
{
    return this->Cpi(celli)*alpha_[celli];
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::kappaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New
    (
        "kappaEff",
        Cp()*(alpha_ + alphat)
    );
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return Cp(patchi)*(alpha_.boundaryField()[patchi] + alphat);
}
// ************************************************************************* //
