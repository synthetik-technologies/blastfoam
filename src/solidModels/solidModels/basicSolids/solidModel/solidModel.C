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
    defineTypeNameAndDebug(solidModel, 0);
    defineRunTimeSelectionTable(solidModel, dictionary);
    defineRunTimeSelectionTable(solidModel, lagrangian);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidModel::checkWedges() const
{
    const fvMesh& mesh = this->mesh();

    label nWedgePatches = 0;
    vector wedgeDirVec = vector::zero;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchI]
            );

            nWedgePatches++;
            wedgeDirVec += cmptMag(wpp.centreNormal());

            // Make sure that solidWedge is used instead of wedge
            if
            (
                DD_.boundaryField()[patchI].type() == "wedge"
             && D_.boundaryField()[patchI].type() == "wedge"
            )
            {
                FatalErrorIn("void Foam::solidModel::checkWedges() const")
                    << "solidWedge should be used on displacement solution "
                    << "field wedge patches as non-orthogonal corrections "
                    << "are important!"
                    << abort(FatalError);
            }
        }
    }

    reduce(nWedgePatches, maxOp<label>());

    if (nWedgePatches)
    {
        if (nWedgePatches != 2)
        {
            FatalErrorIn("void Foam::solidModel::checkWedges() const")
                << "For axisymmetric cases, there should be exactly two wedge "
                << "patches!" << abort(FatalError);
        }

        Info<< nl << "Axisymmetric case: disabling the solution in the "
            << "out-of-plane direction" << endl;

        // We will use const_cast to disable the out-of-lane direction
        Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec /= mag(wedgeDirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                solD[cmpt] = -1;

                wordList dirs(3);
                dirs[0] = "x";
                dirs[1] = "y";
                dirs[2] = "z";
                Info<< "    out-of-plane direction: " << dirs[cmpt] << nl
                    << endl;
            }
            else
            {
                solD[cmpt] = 1;
            }
        }
    }


    // Check all the face normals are in the same direction on the wedge patches
    // This is to avoid the case where a wedge patch is composed of two
    // disconnected regions with one on the front and one on the back
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            // Unit face normals on processor
            const vectorField nf = mesh.boundaryMesh()[patchI].faceNormals();

            if (nf.size() == 0)
            {
                FatalErrorIn("void Foam::solidModel::checkWedges() const")
                    << "There are no faces on the wedge patch "
                    << mesh.boundaryMesh()[patchI].name() << " on this processor:"
                    << nl << "Every processor should have at least one face on "
                    << "each wedge patch"
                    << abort(FatalError);
            }

            // Check that all the wedge face normals point in the same direction

            vector firstFaceNOnMasterProc = vector::zero;

            if (Pstream::master())
            {
                firstFaceNOnMasterProc = nf[0];
            }

            // Sync in parallel so that all processors have the master vector
            reduce(firstFaceNOnMasterProc, sumOp<vector>());

            forAll(nf, faceI)
            {
                if ((nf[faceI] & firstFaceNOnMasterProc) < 0)
                {
                    FatalErrorIn("void Foam::solidModel::checkWedges() const")
                        << "On wedge patch "
                        << mesh.boundaryMesh()[patchI].name()
                        << " there are at "
                        << "least two faces with unit normals in the opposite "
                        << "directions" << nl
                        << "Please check that the wedge patches are correctly "
                        << "defined"
                        << abort(FatalError);
                }
            }
        }
    }
}


const Foam::pointVectorField& Foam::solidModel::pointDorPointDD() const
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


void Foam::solidModel::makeSetCellDisps() const
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


const Foam::setCellDisplacements& Foam::solidModel::setCellDisps() const
{
    if (setCellDispsPtr_.empty())
    {
        makeSetCellDisps();
    }

    return setCellDispsPtr_();
}


// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

Foam::thermalModel& Foam::solidModel::thermal()
{
    return thermal_;
}


Foam::mechanicalModel& Foam::solidModel::mechanical()
{
    return mechanical_;
}


Foam::volScalarField& Foam::solidModel::rho()
{
    return thermal_.rho();
}


void Foam::solidModel::setCellDisps(fvVectorMatrix& DEqn)
{
    if (setCellDisps().cellIDs().size() > 0)
    {
        DEqn.setValues(setCellDisps().cellIDs(), setCellDisps().cellDisps());
    }
}


Foam::wordList Foam::solidModel::pointDBoundaryTypes
(
    const volVectorField& D
) const
{
    wordList bTypes
    (
        D.boundaryField().size(),
        "calculated"
    );
    forAll(D.boundaryField(), patchi)
    {
        if (isA<fixedValueFvPatchVectorField>(D.boundaryField()[patchi]))
        {
            bTypes[patchi] = "fixedValue";
        }
    }
    return bTypes;
}


Foam::dictionary& Foam::solidModel::solidModelDict()
{
    if (this->isDict(type_ + "Coeffs"))
    {
        return this->subDict(type_ + "Coeffs");
    }
    return *this;
}


void Foam::solidModel::displacementFromVelocity
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

void Foam::solidModel::readDict()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModel::solidModel(const word& type, fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "solidProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    ),
    mesh_(mesh),
    type_(type),
    mechanical_(mesh),
    thermal_(mesh, true),
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
    enforceLinear_(false),
    globalPatches_(globalPolyBoundaryMesh::New(mesh))
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


Foam::solidModel::solidModel
(
    const word& type,
    fvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonlinear,
    const bool incremental,
    const bool isSolid
)
:
    IOdictionary
    (
        IOobject
        (
            "solidProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    type_(type),
    mechanical_(mesh, nonlinear, incremental),
    thermal_(mesh, isSolid),
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
    stabilisationPtr_(new momentumStabilisation(solidModelDict())),
    enforceLinear_(false),
    globalPatches_(globalPolyBoundaryMesh::New(mesh))
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

Foam::solidModel::~solidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidModel::initialize()
{
    if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        globalPatches_.setDisplacementField(mesh_.name(), "pointD");
    }
    else
    {
        globalPatches_.setDisplacementField(mesh_.name(), "none");
    }
    globalPatches_.setInverseDisplacement(this->mesh().name(), false);
    globalPatches_.update();
}

const Foam::volScalarField& Foam::solidModel::rho() const
{
    return thermal_.rho();
}

const Foam::thermalModel& Foam::solidModel::thermal() const
{
    return thermal_;
}


const Foam::mechanicalModel& Foam::solidModel::mechanical() const
{
    return mechanical_;
}


void Foam::solidModel::DisRequired(const word& type)
{
    if (!D().headerOk())
    {
        FatalErrorInFunction
            << type << " requires the 'D' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::DDisRequired(const word& type)
{
    if (!DD().headerOk())
    {
        FatalErrorInFunction
            << type << " requires the 'DD' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::updateTotalFields()
{
    thermal().correct();
    mechanical().updateTotalFields();
}


Foam::tmp<Foam::vectorField> Foam::solidModel::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    const symmTensorField& psigma = this->sigma(patch);
    vectorField n(this->nf(patch));

    // Return patch snGrad
    return
        (traction - n*pressure - (n & psigma))/this->impK(patch)
      + (patch.nf() & (this->solutionGradD().boundaryField()[patch.index()]));
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const volScalarField& J)
{
    scalar minJ = min(J).value();
    scalar maxJ = max(J).value();
    if ((minJ < 0.01) || (maxJ > 100))
    {
        DebugInfo<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const surfaceScalarField& J)
{
    scalar minJ = min(J).value();
    scalar maxJ = max(J).value();
    if ((minJ < 0.01) || (maxJ > 100))
    {
        DebugInfo<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


bool Foam::solidModel::read()
{
    if (regIOobject::read())
    {
        readDict();

        return true;
    }
    else
    {
        return false;
    }
}


// bool Foam::solidModel::readIfModified()
// {
//     if (regIOobject::readIfModified())
//     {
//         // Clear current settings except fluxRequired
//         readDict();
//
//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }


bool Foam::solidModel::write(const bool write) const
{
    return true;
}


bool Foam::solidModel::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    return this->write(write);
}


Foam::scalar Foam::solidModel::newDeltaT()
{
    return min
    (
        runTime().deltaTValue(),
        mechanical().newDeltaT()
    );
}

void Foam::solidModel::moveMesh
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


const Foam::dictionary& Foam::solidModel::solidModelDict() const
{
    return this->optionalSubDict(type_ + "Coeffs");
}


// ************************************************************************* //
