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

#include "DSolidModel.H"
#include "volFields.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
#include "RectangularMatrix.H"
#include "solidTractionFvPatchVectorField.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "calculatedPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "globalPolyBoundaryMesh.H"
#include "fvm.H"

#include "fvcGradf.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DSolidModel, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pointVectorField& Foam::DSolidModel::pointDorPointDD() const
{
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        // Updated Lagrangian approaches move the mesh at the end of each
        // time-step so we use the increment of displacement field to calculate
        // the current deformed face zone points
        return pointDD();
    }
    else
    {
        // As linearGeometry and total Lagrangian approaches do not move the
        // mesh, we use the total displacement field to calculate the current
        // deformed face zone points
        return pointD();
    }
}


void Foam::DSolidModel::makeSetCellDisps() const
{
    if (setCellDispsPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set!"
            << abort(FatalError);
    }

    if (this->found("cellDisplacements"))
    {
        setCellDispsPtr_.set
        (
            new setCellDisplacements(mesh(), *this)
        );
    }
    else
    {
        setCellDispsPtr_.set(new setCellDisplacements(mesh()));
    }
}


const Foam::setCellDisplacements& Foam::DSolidModel::setCellDisps() const
{
    if (setCellDispsPtr_.empty())
    {
        makeSetCellDisps();
    }

    return setCellDispsPtr_();
}


// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //


Foam::mechanicalModel& Foam::DSolidModel::mechanical()
{
    return mechanical_;
}


void Foam::DSolidModel::setCellDisps(fvVectorMatrix& DEqn)
{
    if (setCellDisps().cellIDs().size() > 0)
    {
        DEqn.setValues(setCellDisps().cellIDs(), setCellDisps().cellDisps());
    }
}


void Foam::DSolidModel::displacementFromVelocity
(
    volVectorField& disp,
    volVectorField& ddisp
)
{
    // Get the number of stored times
    // We don't care if D or DD is used since both need to be consistent
    label nOld = max(disp.nOldTimes(), ddisp.nOldTimes());

    // Return if no old times are present or we are continuing a simulation
    if (mesh_.time().restart() || !nOld)
    {
        return;
    }

    // If U is found calculate D and DD
    if (!U().headerOk())
    {
        ddisp.primitiveFieldRef() = U()*mesh_.time().deltaT();
        disp.oldTimeRef().primitiveFieldRef() -= ddisp;

        // Call this function on the old times
        displacementFromVelocity
        (
            disp.oldTimeRef(),
            ddisp.oldTimeRef()
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DSolidModel::DSolidModel(const word& type, fvMesh& mesh)
:
    solidModel(type, mesh, true),
    mechanical_(mesh),
    D_
    (
        IOobject
        (
            "D",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    DD_
    (
        IOobject
        (
            "DD",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength/dimTime, vector::zero)
    ),
    pointD_
    (
        IOobject
        (
            "pointD",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimLength, Zero)
    ),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimLength, Zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, Zero)
    ),
    stabilisationPtr_(new momentumStabilisation(solidModelDict())),
    enforceLinear_(false)
{
    globalPatches_.setDisplacementField(mesh_.name(), "none");

    if (!pointD_.headerOk())
    {
        mechanical().volToPoint().interpolate(D_, pointD_);
    }
    if (!pointDD_.headerOk())
    {
        mechanical().volToPoint().interpolate(DD_, pointDD_);
    }
}


Foam::DSolidModel::DSolidModel
(
    const word& type,
    fvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonlinear,
    const bool incremental,
    const bool isSolid
)
:
    solidModel(type, mesh, isSolid),
    mechanical_(mesh, nonlinear, incremental),
    D_
    (
        IOobject
        (
            "D",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    DD_
    (
        IOobject
        (
            "DD",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength/dimTime, vector::zero)
    ),
    pointD_
    (
        IOobject
        (
            "pointD",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimLength, Zero),
        pointDBoundaryTypes(incremental ? DD_ : D_)
    ),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimLength, Zero),
        pointDBoundaryTypes(incremental ? DD_ : D_)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, Zero)
    ),
    stabilisationPtr_(new momentumStabilisation(solidModelDict()))
{
    globalPatches_.setDisplacementField(mesh_.name(), "none");

    if (!pointD_.headerOk())
    {
        mechanical().volToPoint().interpolate(D_, pointD_);
    }
    if (!pointDD_.headerOk())
    {
        mechanical().volToPoint().interpolate(DD_, pointDD_);
    }


    // If the case is axisymmetric, we will disable solving in the out-of-plane
    // direction
    // PC, 12-Nov-18: disabling the 3rd direction slows down convergence a lot
    // in some elastic cases: disabled for now
    //checkWedges();

    if
    (
        solidModelDict().lookupOrDefault
        (
            "initializeDisplacementFromVelocity",
            false
        )
     && !mesh.time().restart()
    )
    {
        typeIOobject<volVectorField> omegaIO
        (
            "omega",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        if (omegaIO.headerOk())
        {
            dimensionedVector xc
            (
                "centreOfRotation",
                dimLength,
                solidModelDict()
            );
            volVectorField omega(omegaIO, mesh_);
            U_ = omega ^ (mesh_.C() - xc);
        }
        else if (solidModelDict().found("omega"))
        {
            dimensionedVector xc
            (
                "centreOfRotation",
                dimLength,
                solidModelDict()
            );
            dimensionedVector omega("omega", inv(dimTime), solidModelDict());
            U_ = omega ^ (mesh_.C() - xc);
        }

        displacementFromVelocity(D_, DD_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DSolidModel::~DSolidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DSolidModel::initialize()
{
    if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        globalPatches_.setDisplacementField(mesh_.name(), pointD_.name());
    }
    else
    {
        globalPatches_.setDisplacementField(mesh_.name(), "none");
    }
    globalPatches_.setInverseDisplacement(this->mesh().name(), false);
    globalPatches_.update();
}


const Foam::mechanicalModel& Foam::DSolidModel::mechanical() const
{
    return mechanical_;
}


void Foam::DSolidModel::updateTotalFields()
{
    thermal().correct();
    mechanical().updateTotalFields();
}


Foam::scalar Foam::DSolidModel::newDeltaT() const
{
    return min
    (
        runTime().deltaTValue(),
        mechanical().newDeltaT()
    );
}

void Foam::DSolidModel::moveMesh
(
    const pointField& oldPoints,
    const volVectorField& DD,
    pointVectorField& pointDD
)
{
    Info<< "Moving the mesh to the deformed configuration" << nl << endl;

    //- Move mesh by interpolating displacement field to vertices

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD, pointDD);

    // Fix, AW/PC, 22-Dec-20,
    // correctBoundaryConditions should not be called as it causes (global?)
    // points to become out of sync. This results in the error "face area does
    // not match neighbour..."
    //pointDD.correctBoundaryConditions();

    vectorField& pointDDI = pointDD.primitiveFieldRef();

    vectorField newPoints(oldPoints);

    // Correct symmetryPlane points
    vector validD((vector(mesh().geometricD()) + vector::one)/2.0);
    vector invalidD(vector::one - validD);
    forAll(mesh().boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            if
            (
                returnReduce(mesh().boundaryMesh()[patchI].size(), sumOp<int>())
            )
            {
                continue;
            }

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());

            const vector i(1, 0, 0);
            const vector j(0, 1, 0);
            const vector k(0, 0, 1);

            if (mag(avgN & i) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].x() = 0;
                }
            }
            else if (mag(avgN & j) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].y() = 0;
                }
            }
            else if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].z() = 0;
                }
            }
        }
        else if (isA<emptyPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            if (!returnReduce(meshPoints.size(), sumOp<int>()))
            {
                continue;
            }

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());

            if (mag(avgN & invalidD) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]] =
                        cmptMultiply(pointDDI[meshPoints[pI]], validD);
                }
            }
        }
    }

    // Note: allPoints will have more points than pointDD if there are
    // globalFaceZones
    forAll(pointDDI, pointI)
    {
        newPoints[pointI] += pointDDI[pointI];
    }

    // Move unused globalFaceZone points
    // Not need anymore as globalFaceZones are not used
    //updateGlobalFaceZoneNewPoints(pointDDI, newPoints);

    const twoDPointCorrector& corrector = twoDPointCorrector::New(mesh());
    corrector.correctPoints(newPoints);
    corrector.correctPoints(pointDDI);
    mesh().movePoints(newPoints);
    mesh().V00();
    // mesh().moving(false);
    const_cast<surfaceScalarField&>(mesh().phi()).writeOpt() =
        IOobject::NO_WRITE;
}


// ************************************************************************* //
