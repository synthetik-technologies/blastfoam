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

#include "basicBlastThermo.H"
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
    defineTypeNameAndDebug(basicBlastThermo, 0);
}

// * * * * * * * * * * * * * * Protected Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicBlastThermo::lookupOrConstruct
(
    const fvMesh& mesh,
    const word& name,
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


Foam::volScalarField& Foam::basicBlastThermo::lookupOrConstructE
(
    const fvMesh& mesh,
    const word& name
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

Foam::basicBlastThermo::basicBlastThermo
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
    integrationSystem
    (
        IOobject::groupName("basicThermo", phaseName),
        p.mesh()
    ),
    thermoDict_(dict),
    master_(master),
    masterName_(masterName),
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
    limit_(dict.lookupOrDefault("limit", true)),
    residualAlpha_("residualAlpha", dimless, 0.0),
    residualRho_("residualRho", dimDensity, 0.0)
{}


Foam::basicBlastThermo::basicBlastThermo
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master,
    const word& masterName
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
    integrationSystem
    (
        IOobject::groupName("basicThermo", phaseName),
        mesh
    ),
    thermoDict_(dict),
    master_(master),
    masterName_(masterName),
    name_(phaseName),
    p_
    (
        lookupOrConstruct
        (
            mesh,
            IOobject::groupName("p", masterName),
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
            IOobject::groupName("rho", phaseName),
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
            IOobject::groupName("T", masterName),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            dimTemperature
        )
    ),
    e_
    (
        lookupOrConstructE
        (
            mesh,
            IOobject::groupName("e", masterName)
        )
    ),
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
    limit_(dict.lookupOrDefault("limit", true)),
    residualAlpha_("residualAlpha", dimless, 0.0),
    residualRho_("residualRho", dimDensity, 0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicBlastThermo::~basicBlastThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicBlastThermo::correct()
{
    if (master_)
    {
        this->T_ = this->calcT();
        this->T_.correctBoundaryConditions();
    }
}

Foam::word Foam::basicBlastThermo::readThermoType(const dictionary& dict)
{
    return word
    (
        word(dict.lookup("transport")) + '<'
      + word(dict.lookup("thermo")) + '<'
      + word(dict.lookup("equationOfState"))  + '<'
      + word("specie") + ">>>"
    );
}

Foam::wordList Foam::basicBlastThermo::splitThermoName
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


Foam::wordList Foam::basicBlastThermo::splitThermoName
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
Foam::basicBlastThermo::eBoundaryTypes(const volScalarField& T)
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
Foam::basicBlastThermo::eBoundaryBaseTypes(const volScalarField& T)
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



void Foam::basicBlastThermo::eBoundaryCorrection()
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


Foam::tmp<Foam::volScalarField> Foam::basicBlastThermo::Y(const word& name) const
{
    return volScalarField::New
    (
        IOobject::groupName(name, this->name()),
        e_.mesh(),
        dimless
    );
}


Foam::tmp<Foam::volScalarField> Foam::basicBlastThermo::Y(const label i) const
{
    return volScalarField::New
    (
        IOobject::groupName("Yi" + Foam::name(i), this->name()),
        e_.mesh(),
        dimless
    );
}


void Foam::basicBlastThermo::addDelta
(
    const word& name,
    const volScalarField& delta
)
{}


Foam::tmp<Foam::volScalarField>
Foam::basicBlastThermo::kappa() const
{
    tmp<Foam::volScalarField> kappa(Cp()*this->alpha_);
    kappa.ref().rename("kappa");
    return kappa;
}


Foam::tmp<Foam::scalarField> Foam::basicBlastThermo::kappa
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


Foam::scalar Foam::basicBlastThermo::kappai
(
    const label celli
) const
{
    return this->Cpi(celli)*this->alpha_[celli];
}


Foam::tmp<Foam::volScalarField>
Foam::basicBlastThermo::alphahe() const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCv()*this->alpha_);
    alphaEff.ref().rename("alphahe");
    return alphaEff;
}


Foam::tmp<Foam::scalarField>
Foam::basicBlastThermo::alphahe(const label patchi) const
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
Foam::basicBlastThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> kappaEff(Cp()*(this->alpha_ + alphat));
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


Foam::tmp<Foam::scalarField>
Foam::basicBlastThermo::kappaEff
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
Foam::basicBlastThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCv()*(this->alpha_ + alphat));
    alphaEff.ref().rename("alphaEff");
    return alphaEff;
}


Foam::tmp<Foam::scalarField>
Foam::basicBlastThermo::alphaEff
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
