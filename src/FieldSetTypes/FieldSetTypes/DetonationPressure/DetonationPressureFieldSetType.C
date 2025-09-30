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

#include "DetonationPressureFieldSetType.H"
#include "hashedWordList.H"
#include "fluidBlastThermo.H"
#include "FieldSetTypesFwd.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace FieldSetTypes
{
    defineTypeNameAndDebug(DetonationPressure, 0);
    addToRunTimeSelectionTable(scalarVolFieldSetType, DetonationPressure, dictionary);

}
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::FieldSetTypes::DetonationPressure::DetonationPressure
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const labelList& selectedIndices,
    Istream& is,
    const bool write
)
:
    VolFieldSetType<scalar>
    (
        mesh,
        dict,
        fieldName,
        selectedIndices,
        is,
        write
    ),
    phaseName_(is),
    e0_(readScalar(is)),
    pDetPtr_(nullptr)
{
    if (fieldName != "p")
    {
        FatalErrorInFunction
            << typeName << " should only be used to set pressure" << endl
            << abort(FatalError);
    }
    if (is.good())
    {
        e0_ /= readScalar(is);
    }

    if
    (
        lookupOrRead<volScalarField>(IOobject::groupName("alpha", phaseName_))
     && this->good_
    )
    {
        // If rho is already found, check out to remove conflicts
        if (mesh.foundObject<volScalarField>("rho"))
        {
            mesh.lookupObjectRef<volScalarField>("rho").checkOut();
        }

        fluidBlastThermo& thermo = lookupOrConstructThermo(mesh, word::null);
        // Make sure this rho field is the one in the database
        thermo.rhoRef().checkIn();

        // Do not write energy
        thermo.he().writeOpt() = IOobject::NO_WRITE;
        pDetPtr_ = this->lookupOrConstruct
        (
            IOobject::groupName("pDet", phaseName_),
            thermo.p()
        );
        volScalarField& pDet = *pDetPtr_;
        pDet = thermo.p();

        const volScalarField& alpha =
            *lookupOrRead<volScalarField>
            (
                IOobject::groupName("alpha", phaseName_)
            );
        volScalarField& he = thermo.he();
        forAll(alpha, celli)
        {
            const scalar heOld = he[celli];
            he[celli] += e0_*alpha[celli];
            pDet[celli] = thermo.cellpRhoT(celli);
            he[celli] = heOld;
        }
        this->setField();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FieldSetTypes::DetonationPressure::~DetonationPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word
Foam::FieldSetTypes::DetonationPressure::readPhaseName
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


Foam::dictionary
Foam::FieldSetTypes::DetonationPressure::thermoDict
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

    dictionary phaseDict
    (
        phaseProperties.subDict
        (
            phaseName == word::null
          ? "mixture"
          : phaseName
        )
    );
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


Foam::fluidBlastThermo& Foam::FieldSetTypes::DetonationPressure::lookupOrConstructThermo
(
    const fvMesh& mesh,
    const word& phaseName
) const
{
    word thermoName(IOobject::groupName(physicalProperties::typeName, phaseName));

    if (mesh.foundObject<fluidBlastThermo>(thermoName))
    {
        return mesh.lookupObjectRef<fluidBlastThermo>(thermoName);
    }

    fluidBlastThermo* thermoPtr = nullptr;
    dictionary dict(thermoDict(mesh, phaseName));
    word stateType(dict.lookupOrDefault<word>("stateType", "fluid"));


    if (stateType == "fluid")
    {
        thermoPtr =
            fluidBlastThermo::New
            (
                mesh,
                dict,
                word::null,
                phaseName
            ).ptr();
    }
    else
    {
        FatalErrorInFunction
            << "Only fluids are allowed" << nl
            << abort(FatalError);
    }
    thermoPtr->store(thermoPtr);

    return mesh.lookupObjectRef<fluidBlastThermo>(thermoName);
}


void Foam::FieldSetTypes::DetonationPressure::setGeoField
(
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<scalar>& f,
    const label patchi
)
{
    const volScalarField& pDet = *pDetPtr_;
    if (patchi < 0)
    {
        forAll(indices, i)
        {
            f[i] = pDet[indices[i]];
        }
    }
    else
    {
        forAll(indices, i)
        {
            f[i] = pDet.boundaryField()[patchi][indices[i]];
        }
    }
}


// ************************************************************************* //
