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
        mesh_,
        dimensionedSymmTensor(dimViscosity*dimVelocity/dimLength, Zero)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectListSolver::immersedBoundaryObjectListSolver
(
    const fvMesh& mesh
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
    ibmDict_
    (
        IOobject
        (
            "immersedBoundaryProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE // Handles internally
        )
    ),
    objects_(0),
    test_(false),
    curTimeIndex_(-1),
    cellTypes_
    (
        IOobject
        (
            "cellType",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        0.0
    ),
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

    const dictionary& objectDict(ibmDict_.subDict("objects"));
    wordList objects(objectDict.toc());
    objects_.resize(objects.size());
    bool moving = false;

    dictionary stateDict;
    if (stateDictIO.typeHeaderOk<IOdictionary>(true))
    {
        stateDict = IOdictionary(stateDictIO);
    }
    else
    {
        stateDict = objectDict;
    }
Info<<objects<<endl;
    forAll(objects, i)
    {
        Info<<i<<endl;
        objects_.set
        (
            i,
            objects[i],
            immersedBoundaryObject::New
            (
                mesh,
                objectDict.subDict(objects[i]),
                stateDict.subDict(objects[i])
            ).ptr()
        );
        objects_[i].initialize();
        objects_[i].setInternal(cellTypes_, 2.0);
        objects_[i].setBoundary(cellTypes_, 1.0);
        thermalForcingNeeded_ =
            thermalForcingNeeded_ || objects_[i].temperatureDependent();
        moving = moving || objects_[i].moving();
        Info<<i<<endl;
        Info<<objects_[i].name()<<endl;
    }
    if (moving)
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
    if (objects_.size() > 1)
    {
        objectIDPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "objectID",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless, 0.0)
            )
        );
    }
    Info<<objects_[0].name()<<endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectListSolver::~immersedBoundaryObjectListSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::PtrListDictionary<Foam::immersedBoundaryObject>&
Foam::immersedBoundaryObjectListSolver::objects() const
{
    Info<<objects_[0].name()<<endl;
    return objects_;
}


void Foam::immersedBoundaryObjectListSolver::solve()
{
    cellTypes_ = Zero;
    if (objectIDPtr_.valid())
    {
        objectIDPtr_() = -1.0;
    }

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
        objects_[i].setInternal(cellTypes_, 2.0);
        objects_[i].setBoundary(cellTypes_, 1.0);
        if (objectIDPtr_.valid())
        {
            objects_[i].setInternal(objectIDPtr_(), scalar(i + 1));
            objects_[i].setBoundary(objectIDPtr_(), scalar(i + 1));
        }

        // Zero forces
        objects_[i].force() = Zero;
        objects_[i].forceExt() = Zero;
        objects_[i].momentExt() = Zero;

        objects_[i].force() +=
            objects_[i].Sf()*objects_[i].patchExternalField(p)
          + (objects_[i].Sf() & objects_[i].patchExternalField(sigma));

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

        objects_[i].movePoints();

        if (objects_[i].report())
        {
            objects_[i].status();
        }
    }

//     correctBoundaryConditions();
}


void Foam::immersedBoundaryObjectListSolver::correctBoundaryConditions()
{
    forAll(objects_, i)
    {
        objects_[i].setValues();
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
    Info<<"write"<<endl;
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
