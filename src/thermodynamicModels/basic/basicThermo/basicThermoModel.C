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

#include "basicThermoModel.H"
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
    defineTypeNameAndDebug(basicThermoModel, 0);
}

// * * * * * * * * * * * * * * Protected Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicThermoModel::lookupOrConstruct
(
    const fvMesh& mesh,
    const char* name,
    const IOobject::readOption rOpt,
    const IOobject::writeOption wOpt,
    const dimensionSet& dims
) const
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
    {
        volScalarField* fPtr;

        if (rOpt == IOobject::MUST_READ)
        {
            fPtr =
                new volScalarField
                (
                    IOobject
                    (
                        name,
                        mesh.time().timeName(),
                        mesh,
                        rOpt,
                        wOpt
                    ),
                    mesh
                );
        }
        else
        {
            fPtr =
                new volScalarField
                (
                    IOobject
                    (
                        name,
                        mesh.time().timeName(),
                        mesh,
                        rOpt,
                        wOpt
                    ),
                    mesh,
                    dimensionedScalar("0", dims, Zero)
                );
        }

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volScalarField>(name);
}


Foam::volScalarField& Foam::basicThermoModel::lookupOrConstructE
(
    const fvMesh& mesh,
    const char* name
) const
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
    {
        volScalarField* fPtr
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimEnergy/dimMass, -1.0),
                eBoundaryTypes(T_),
                eBoundaryBaseTypes(T_)
            )
        );
        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volScalarField>(name);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicThermoModel::basicThermoModel
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
            IOobject::groupName("basicThermo", phaseName),
            p.mesh().time().timeName(),
            p.mesh()
        )
    ),
    master_(master),
    name_(phaseName),
    p_(p),
    rho_(rho),
    T_(T),
    e_(e),
    alpha_
    (
        IOobject
        (
            IOobject::groupName("thermo:alpha", name_),
            p_.mesh().time().timeName(),
            p_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    limit_(dict.lookupOrDefault("limit", true))
{}


Foam::basicThermoModel::basicThermoModel
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName("basicThermo", phaseName),
            mesh.time().timeName(),
            mesh
        )
    ),
    master_(master),
    name_(phaseName),
    p_
    (
        lookupOrConstruct
        (
            mesh,
            "p",
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            dimPressure
        )
    ),
    rho_
    (
        lookupOrConstruct
        (
            mesh,
            "rho",
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            dimDensity
        )
    ),
    T_
    (
        lookupOrConstruct
        (
            mesh,
            "T",
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            dimTemperature
        )
    ),
    e_(lookupOrConstructE(mesh, "e")),
    alpha_
    (
        IOobject
        (
            IOobject::groupName("thermo:alpha", name_),
            p_.mesh().time().timeName(),
            p_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    limit_(dict.lookupOrDefault("limit", true))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermoModel::~basicThermoModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicThermoModel::correct()
{
    if (master_)
    {
        T_ = calcT();
    }
}

Foam::word Foam::basicThermoModel::readThermoType(const dictionary& dict)
{
    return word
    (
        word(dict.lookup("transport")) + '<'
      + word(dict.lookup("thermo")) + '<'
      + word(dict.lookup("equationOfState"))  + '<'
      + word("specie") + ">>>"
    );
}


Foam::wordList Foam::basicThermoModel::splitThermoName
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
Foam::basicThermoModel::eBoundaryTypes(const volScalarField& T)
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
Foam::basicThermoModel::eBoundaryBaseTypes(const volScalarField& T)
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



void Foam::basicThermoModel::eBoundaryCorrection()
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


Foam::tmp<Foam::volScalarField>
Foam::basicThermoModel::kappa() const
{
    tmp<Foam::volScalarField> kappa(Cp()*this->alpha_);
    kappa.ref().rename("kappa");
    return kappa;
}


Foam::tmp<Foam::scalarField> Foam::basicThermoModel::kappa
(
    const label patchi
) const
{
    return
        Cp
        (
            this->rho_.boundaryField()[patchi],
            this->e_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}


Foam::scalar Foam::basicThermoModel::kappai
(
    const label celli
) const
{
    return this->Cpi(celli)*this->alpha_[celli];
}


Foam::tmp<Foam::volScalarField>
Foam::basicThermoModel::alphahe() const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCv()*this->alpha_);
    alphaEff.ref().rename("alphahe");
    return alphaEff;
}


Foam::tmp<Foam::scalarField>
Foam::basicThermoModel::alphahe(const label patchi) const
{
    return
        this->CpByCv
        (
            this->rho_.boundaryField()[patchi],
            this->e_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField>
Foam::basicThermoModel::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> kappaEff(Cp()*(this->alpha_ + alphat));
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


Foam::tmp<Foam::scalarField>
Foam::basicThermoModel::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->Cp
        (
            this->rho_.boundaryField()[patchi],
            this->e_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )
       *(
           this->alpha_.boundaryField()[patchi]
         + alphat
        );
}


Foam::tmp<Foam::volScalarField>
Foam::basicThermoModel::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCv()*(this->alpha_ + alphat));
    alphaEff.ref().rename("alphaEff");
    return alphaEff;
}


Foam::tmp<Foam::scalarField>
Foam::basicThermoModel::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->CpByCv
        (
            this->rho_.boundaryField()[patchi],
            this->e_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )
       *(
            this->alpha_.boundaryField()[patchi]
          + alphat
        );
}

// ************************************************************************* //
