/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "solidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidModel, 0);
    defineRunTimeSelectionTable(solidModel, dictionary);
    defineRunTimeSelectionTable(solidModel, lagrangian);
}


Foam::tmp<Foam::volScalarField> Foam::solidModel::uniformOrRead
(
    const fvMesh& mesh,
    const dictionary& dict,
    dimensionedScalar val
)
{
    const word entry = val.name();
    const dictionary& entryDict = dict.subDict(entry);
    const word type = entryDict.lookup("type");

    IOobject fieldIO
    (
        entry,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    if (type == "uniform")
    {
        val.value() = entryDict.lookup<scalar>("value");
        return tmp<volScalarField>
        (
            new volScalarField
            (
                fieldIO,
                mesh,
                val
            )
        );
    }
    else if (type == "field")
    {
        fieldIO.readOpt() = IOobject::MUST_READ;
        return tmp<volScalarField>
        (
            new volScalarField
            (
                fieldIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for " << entry
            << abort(FatalError);
    }

    return tmp<volScalarField>();
}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

void Foam::solidModel::DisRequired()
{
    IOobject DIO
    (
        "D",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!DIO.typeHeaderOk<volVectorField>(true))
    {
        FatalErrorInFunction
            << "This solidModel requires the 'D' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::DDisRequired()
{
    IOobject DDIO
    (
        "DD",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!DDIO.typeHeaderOk<volVectorField>(true))
    {
        FatalErrorInFunction
            << "This solidModel requires the 'DD' field to be specified!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModel::solidModel
(
    const word& type,
    fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "solidProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    type_(type),
    thermoPtr_(),
    D_
    (
        IOobject
        (
            "D",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    pointD_
    (
        IOobject
        (
            "pointD",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedVector("0", dimLength, Zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    relTol_
    (
        subDict(type + "Coeffs").lookupOrDefault<scalar>
        (
            "relTol",
            1e-04
        )
    ),
    tolerance_
    (
        subDict(type + "Coeffs").lookupOrDefault<scalar>
        (
            "solutionTolerance",
            1e-06
        )
    ),
    nCorr_(subDict(type + "Coeffs").lookupOrDefault<label>("nCorrectors", 10000)),
    fvModels_(fvModels::New(mesh)),
    fvConstraints_(fvConstraints::New(mesh))
{
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    thermoPtr_.set(solidBlastThermo::New(mesh, thermophysicalProperties).ptr());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidModel::~solidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
