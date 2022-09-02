/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "MassIntegratedFieldSetType.H"
#include "hashedWordList.H"
#include "fluidBlastThermo.H"
#include "solidBlastThermo.H"
#include "calcAngleFraction.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
Foam::FieldSetTypes::MassIntegrated<Type, FSType>::MassIntegrated
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const labelList& selectedIndices,
    Istream& is,
    const bool write
)
:
    FSType<Type>
    (
        mesh,
        dict,
        fieldName,
        selectedIndices,
        is,
        write
    ),
    phaseName_(readPhaseName(is, fieldName)),
    value_(pTraits<Type>(is)),
    thermo_(lookupOrConstructThermo(mesh, phaseName_))
{
    if (this->good_)
    {
        const volScalarField& rho(thermo_.rho());
        scalar mass(0.0);
        forAll(selectedIndices, i)
        {
            label celli = selectedIndices[i];
            mass += rho[celli]*mesh.V()[celli];
        }
        mass /= calcAngleFraction(mesh);

        this->value_ = this->value_/mass;

        this->setField();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
Foam::FieldSetTypes::MassIntegrated<Type, FSType>::~MassIntegrated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
Foam::word
Foam::FieldSetTypes::MassIntegrated<Type, FSType>::readPhaseName
(
    Istream& is,
    const word& fieldName
) const
{
    word phaseName(IOobject::group(fieldName));
    token t(is);
    if (t.isWord())
    {
        phaseName = t.wordToken();
    }
    else
    {
        is.putBack(t);
    }
    return phaseName;
}


template<class Type, template<class> class FSType>
Foam::dictionary
Foam::FieldSetTypes::MassIntegrated<Type, FSType>::thermoDict
(
    const fvMesh& mesh,
    const word& phaseName
) const
{
    IOdictionary phaseProperties
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        )
    );

    hashedWordList phases
    (
        phaseProperties.lookupOrDefault
        (
            "phases",
            hashedWordList({"mixture"})
        )
    );

    dictionary phaseDict(phaseProperties.subDict(phaseName));
    word stateType;
    if (phaseDict.found("thermoType"))
    {
        stateType =
            phaseDict.subDict("thermoType").lookup<word>("equationOfState")
         == "rhoConst"
          ? "solid"
          : "fluid";
    }
    else if (phaseDict.found("products"))
    {
        stateType =
            phaseDict.subDict("products").subDict("thermoType").lookup<word>
            (
                "equationOfState"
            ) == "rhoConst"
          ? "solid"
          : "fluid";
    }
    else
    {
        FatalErrorInFunction
            << "Could not determine state type" << nl
            << "neither thermoType or products/thermoDict was found in" << nl
            << phaseDict
            << abort(FatalError);
    }
    phaseDict.set("stateType", stateType);
    return phaseDict;
}


template<class Type, template<class> class FSType>
const Foam::blastThermo&
Foam::FieldSetTypes::MassIntegrated<Type, FSType>::lookupOrConstructThermo
(
    const fvMesh& mesh,
    const word& phaseName
) const
{
    word thermoName(IOobject::groupName(basicThermo::dictName, phaseName));

    if (mesh.foundObject<blastThermo>(thermoName))
    {
        return mesh.lookupObject<blastThermo>(thermoName);
    }

    blastThermo* thermoPtr = nullptr;
    dictionary dict(thermoDict(mesh, phaseName));
    word stateType(dict.lookup<word>("stateType"));

    if
    (
        (stateType == "propellant" || stateType == "solid")
     && !mesh.foundObject<volScalarField>
        (
            IOobject::groupName("p", phaseName)
        )
    )
    {
        volScalarField* psPtr =
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("p", phaseName),
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimPressure, 0.0)
            );
        psPtr->store(psPtr);
    }

    if (stateType == "fluid")
    {
        thermoPtr =
            fluidBlastThermo::New
            (
                1,
                mesh,
                dict,
                phaseName
            ).ptr();
    }
    else if (stateType == "solid")
    {
        thermoPtr =
            solidBlastThermo::New
            (
                mesh,
                dict,
                phaseName
            ).ptr();
    }
    else
    {
        FatalErrorInFunction
            << "unknown state " << stateType << nl
            << abort(FatalError);
    }
    thermoPtr->store(thermoPtr);

    return mesh.lookupObject<blastThermo>(thermoName);
}


template<class Type, template<class> class FSType>
void
Foam::FieldSetTypes::MassIntegrated<Type, FSType>::getInternalField
(
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{
    f = value_;
}


template<class Type, template<class> class FSType>
void
Foam::FieldSetTypes::MassIntegrated<Type, FSType>::getBoundaryField
(
    const label patchi,
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{
    f = value_;
}

// ************************************************************************* //
