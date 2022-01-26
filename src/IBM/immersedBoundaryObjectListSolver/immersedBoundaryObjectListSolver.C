/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "immersedBoundaryObjectListSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "uniformDimensionedFields.H"
#include "dimensionedTypes.H"
#include "mathematicalConstants.H"
#include "kinematicMomentumTransportModel.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidThermo.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryObjectListSolver, 0);
}

Foam::tmp<Foam::volSymmTensorField>
Foam::immersedBoundaryObjectListSolver::devTau(const word& phaseName) const
{
    if (mesh_.foundObject<volSymmTensorField>("devTau"))
    {
        return mesh_.lookupObject<volSymmTensorField>("devTau");
    }

    typedef incompressible::momentumTransportModel icoModel;
    typedef compressible::momentumTransportModel cmpModel;

    const word& modelName = momentumTransportModel::typeName;
    if (mesh_.foundObject<icoModel>(modelName))
    {
        const incompressible::momentumTransportModel& model =
            mesh_.lookupObject<icoModel>(modelName);

        return mesh_.lookupObject<volScalarField>("rho")*model.devSigma();
    }
    else if (mesh_.foundObject<cmpModel>(modelName))
    {
        const cmpModel& model =
            mesh_.lookupObject<cmpModel>(modelName);

        return model.devTau();
    }
    else if (mesh_.foundObject<dictionary>("transportProperties"))
    {
        // Legacy support for icoFoam

        const dictionary& transportProperties =
             mesh_.lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        const volVectorField& U =
            mesh_.lookupObject<volVectorField>("U");

        return
          - mesh_.lookupObject<volScalarField>("rho")
           *nu*dev(twoSymm(fvc::grad(U)));
    }
    else if (debug)
    {
        WarningInFunction
            << "No valid model for viscous stress calculation" << endl;
    }

    return volSymmTensorField::New
    (
        "devTau",
        dynamicCast<const fvMesh&>(mesh_),
        dimensionedSymmTensor(dimViscosity*dimVelocity/dimLength, Zero)
    );
}


Foam::volScalarField& Foam::immersedBoundaryObjectListSolver::cellTypes()
{
    if (!cellTypesPtr_.valid())
    {
        cellTypesPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "cellTypes",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                dynamicCast<const fvMesh&>(mesh_),
                0.0
            )
        );
    }
    return cellTypesPtr_();
}


Foam::volScalarField& Foam::immersedBoundaryObjectListSolver::objectIDs()
{
    if (!objectIDPtr_.valid())
    {
        objectIDPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "objectIDs",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                dynamicCast<const fvMesh&>(mesh_),
                -1
            )
        );
    }
    return objectIDPtr_();
}


const Foam::dictionary&
Foam::immersedBoundaryObjectListSolver::ibmProperties
(
    const polyMesh& mesh
) const
{
    const Time& runTime = mesh.time();
    if (!runTime.foundObject<IOdictionary>("immersedBoundaryProperties"))
    {
        IOdictionary* ibmPropertiesPtr =
            new IOdictionary
            (
                IOobject
                (
                    "immersedBoundaryProperties",
                    mesh.time().constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE // Handles internally
                )
            );
        ibmPropertiesPtr->store(ibmPropertiesPtr);
    }

    return
        runTime.lookupObject<IOdictionary>("immersedBoundaryProperties");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectListSolver::immersedBoundaryObjectListSolver
(
    const polyMesh& mesh
)
:
    ImmersedBoundaryObjectListSolverObject
    (
        mesh,
        IOobject
        (
            typeName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    mesh_(mesh),
    ibmDict_(ibmProperties(mesh)),
    objects_(0),
    test_(false),
    curTimeIndex_(-1),
    cellTypesPtr_(),
    objectIDPtr_(),
    thermalForcingNeeded_(false),
    collision_()
{
    IOobject stateDictIO
    (
        "immersedBoundaryObjectMotionState",
        mesh.time().timeName(),
        "uniform",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );
    if (stateDictIO.typeHeaderOk<IOdictionary>(true))
    {
        stateDict_ = IOdictionary(stateDictIO);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectListSolver::~immersedBoundaryObjectListSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::immersedBoundaryObject&
Foam::immersedBoundaryObjectListSolver::addObject(const polyPatch& patch)
{
    const word name = patch.name();

    // If the object in the time object registry, using the same name
    // is not the same data as the owner mesh, look up that mesh,
    // and return the given object.
    // This removes duplicate objects when temporary meshes are created
    if (mesh_.time().foundObject<polyMesh>(mesh_.name()))
    {
        if (&mesh_ != &mesh_.time().lookupObject<polyMesh>(mesh_.name()))
        {
            const polyMesh& mesh =
                mesh_.time().lookupObject<polyMesh>(mesh_.name());
            return
                immersedBoundaryObjectListSolver::New(mesh).objects()[name];
        }
    }

    // Return and existing object
    if (objects_.found(name))
    {
        return objects_[name];
    }

    // No existing object was found so create one
    const dictionary& objectDict
    (
        ibmDict_.subDict("objects").subDict(name)
    );

    if (stateDict_.empty())
    {
        stateDict_ = ibmDict_.subDict("objects");
    }

    label i = objects_.size();
    objects_.resize(i + 1);

    objects_.set
    (
        i,
        name,
        immersedBoundaryObject::New
        (
            patch,
            objectDict,
            stateDict_.subDict(name)
        ).ptr()
    );

    thermalForcingNeeded_ =
        thermalForcingNeeded_ || objects_[i].temperatureDependent();

    bool moving = false;
    forAll(objects_, j)
    {
        if (objects_[j].moving())
        {
            moving = true;
            break;
        }
    }


    if (moving && !collision_.valid() && objects_.size() > 1)
    {
        collision_.set
        (
            collisionModel::New
            (
                mesh_,
                objects_,
                ibmDict_.subDict("collisionModel")
            ).ptr()
        );
    }
    return objects_[i];
}



void Foam::immersedBoundaryObjectListSolver::solve()
{
    const Time& t = mesh_.time();

    // Store the motion state at the beginning of the time-stepbool
    if (curTimeIndex_ != t.timeIndex())
    {
        forAll(objects_, i)
        {
            objects_[i].newTime();
        }
        curTimeIndex_ = t.timeIndex();
    }

    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    volSymmTensorField sigma(devTau(word::null));


    forAll(objects_, i)
    {
        // Zero forces
        objects_[i].forceExt() = Zero;
        objects_[i].momentExt() = Zero;

        objects_[i].force() =
            objects_[i].Sf()*objects_[i].patchExternalField(p)
          + (
                objects_[i].Sf()
              & objects_[i].patchExternalField(sigma)
            );

        objects_[i].update();

    }

    if (collision_.valid())
    {
        collision_->collide();
    }

    forAll(objects_, i)
    {
        objects_[i].update();
    }

    if (collision_.valid())
    {
        collision_->updateExternalForce();
    }

    forAll(objects_, i)
    {
        objects_[i].update();

        objects_[i].solve();

        if (Pstream::master())
        {
            objects_[i].status();
        }
    }
    setCellTypes();
    setObjectIDs();
}


void Foam::immersedBoundaryObjectListSolver::setCellTypes()
{
    volScalarField& cellT(cellTypes());
    cellT == Zero;
    forAll(objects_, i)
    {
        objects_[i].setInternal(cellT, 2.0, true);
        objects_[i].setBoundary(cellT, 1.0);
    }
}


void Foam::immersedBoundaryObjectListSolver::setObjectIDs()
{
    if (objects_.size() < 2)
    {
        return;
    }
    volScalarField& objectID(objectIDs());
    objectID == -1;
    forAll(objects_, i)
    {
        objects_[i].setInternal(objectID, scalar(i + 1), true);
        objects_[i].setBoundary(objectID, scalar(i + 1));
    }
}


bool Foam::immersedBoundaryObjectListSolver::movePoints()
{
    forAll(objects_, i)
    {
        objects_[i].clearOut();
    }
    return true;
}


void Foam::immersedBoundaryObjectListSolver::updateMesh(const mapPolyMesh&)
{
    forAll(objects_, i)
    {
        objects_[i].clearOut();
    }
}


Foam::scalar Foam::immersedBoundaryObjectListSolver::CoNum() const
{
    scalar coNum = 0.0;
    forAll(objects_, i)
    {
        coNum = max(coNum, objects_[i].CoNum());
    }

    return coNum;
}


Foam::scalar Foam::immersedBoundaryObjectListSolver::maxCoNum() const
{
    scalar maxcoNum = great;
    forAll(objects_, i)
    {
        maxcoNum = min(maxcoNum, objects_[i].maxCoNum());
    }

    return maxcoNum;
}


bool Foam::immersedBoundaryObjectListSolver::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    IOdictionary dict
    (
        IOobject
        (
            "immersedBoundaryObjectMotionState",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(objects_, i)
    {
        dictionary subDict;
        objects_[i].write(subDict);
        refCast<dictionary>(dict).add(objects_[i].name(), subDict);
    }

    return
        dict.regIOobject::writeObject
        (
            fmt,
            ver,
            cmp,
            write
        );
}

bool Foam::immersedBoundaryObjectListSolver::write(const bool write) const
{
    IOdictionary dict
    (
        IOobject
        (
            "immersedBoundaryObjectMotionState",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(objects_, i)
    {
        dictionary subDict;
        objects_[i].write(subDict);
        refCast<dictionary>(dict).add(objects_[i].name(), subDict);
    }

    return
        dict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            mesh_.time().writeCompression(),
            true
        );
}


// ************************************************************************* //
