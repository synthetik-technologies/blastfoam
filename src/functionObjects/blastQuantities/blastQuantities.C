/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
20-06-2020 Synthetik Applied Technologies: | Time of arrival implementation
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "blastQuantities.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(blastQuantities, 0);
    addToRunTimeSelectionTable(functionObject, blastQuantities, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::volScalarField&
Foam::functionObjects::blastQuantities::lookupOrCreate
(
    const word& name,
    const dimensionSet& dims
) const
{
    Log << "    Reading/initialising field " << name << endl;

    if (obr_.foundObject<volScalarField>(name))
    {
        return obr_.lookupObjectRef<volScalarField>(name);
    }

    // Store on registry
    volScalarField* fieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                name,
                obr_.time().timeName(),
                obr_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("0", dims, 0.0),
            "zeroGradient"
        )
    );
    fieldPtr->store(fieldPtr);

    return *fieldPtr;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blastQuantities::blastQuantities
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    pRef_("pRef", dimPressure, dict),
    overpressureName_
    (
        IOobject::groupName
        (
            "overpressure",
            IOobject::group(pName_)
        )
    ),
    impulse_
    (
        IOobject
        (
            IOobject::groupName
            (
                "impulse",
                IOobject::group(pName_)
            ),
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT
        ),
        mesh_,
        dimensionedScalar("0", dimPressure*dimTime, 0.0)
    ),
    timeOfArrival_
    (
        IOobject
        (
            IOobject::groupName
            (
                "timeOfArrival",
                IOobject::group(pName_)
            ),
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT
        ),
        mesh_,
        runTime.startTime()
    ),
    positivePhaseDuration_
    (
        IOobject
        (
            IOobject::groupName
            (
                "positivePhaseDuration",
                IOobject::group(pName_)
            ),
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT
        ),
        mesh_,
        dimensionedScalar("0", dimTime, 0.0)
    ),
    negativePhaseDuration_
    (
        IOobject
        (
            IOobject::groupName
            (
                "negativePhaseDuration",
                IOobject::group(pName_)
            ),
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT
        ),
        mesh_,
        dimensionedScalar("0", dimTime, 0.0)
    ),
    stage_
    (
        IOobject
        (
            IOobject::groupName
            (
                "blastQuantities:stage",
                IOobject::group(pName_)
            ),
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT
        ),
        mesh_,
        FriedlanderStage::PRE
    ),
    pMax_
    (
        lookupOrCreate
        (
            IOobject::groupName
            (
                IOobject::member(pName_) + "Max",
                IOobject::group(pName_)
            ),
            dimPressure
        )
    )
{
    if (!dict.lookupOrDefault("executeAtStart", false))
    {
        executeAtStart_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::blastQuantities::~blastQuantities()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::blastQuantities::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::blastQuantities::execute()
{
    const fvMesh& mesh = this->mesh_;
    const volScalarField& p(mesh.lookupObject<volScalarField>(pName_));

    // Store overpressure
    store(overpressureName_, p - pRef_);

    // Update impulse
    impulse_ = impulse_.oldTime() + (p - pRef_)*obr_.time().deltaT();

    forAll(pMax_, celli)
    {
        update
        (
            p[celli],
            pMax_[celli],
            stage_[celli],
            timeOfArrival_[celli],
            positivePhaseDuration_[celli],
            negativePhaseDuration_[celli]
        );
    }

    const volScalarField::Boundary& pp = p.boundaryField();
    volScalarField::Boundary& ppMax = pMax_.boundaryFieldRef();
    volScalarField::Boundary& pstage = stage_.boundaryFieldRef();
    volScalarField::Boundary& ptimeOfArrival =
        timeOfArrival_.boundaryFieldRef();
    volScalarField::Boundary& ppositivePhaseDuration =
        positivePhaseDuration_.boundaryFieldRef();
    volScalarField::Boundary& pnegativePhaseDuration =
        negativePhaseDuration_.boundaryFieldRef();
    forAll(pp, patchi)
    {
        forAll(pp[patchi], facei)
        {
            update
            (
                pp[patchi][facei],
                ppMax[patchi][facei],
                pstage[patchi][facei],
                ptimeOfArrival[patchi][facei],
                ppositivePhaseDuration[patchi][facei],
                pnegativePhaseDuration[patchi][facei]
            );
        }
    }

    return true;
}


bool Foam::functionObjects::blastQuantities::write()
{
    return
        writeObject(overpressureName_)
     && pMax_.write()
     && impulse_.write()
     && timeOfArrival_.write()
     && stage_.write()
     && positivePhaseDuration_.write()
     && negativePhaseDuration_.write();
}


// ************************************************************************* //
