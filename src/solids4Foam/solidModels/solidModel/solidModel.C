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
#include "primitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "calculatedPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "fixedValueFvPatchFields.H"

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
        FatalErrorIn(type() + "::makeSetCellDisps() const")
            << "pointer already set!" << abort(FatalError);
    }

    if (solidModelDict().found("cellDisplacements"))
    {
        setCellDispsPtr_.set
        (
            new setCellDisplacements
            (
                mesh(), solidModelDict().subDict("cellDisplacements")
            )
        );
    }
    else
    {
        dictionary dict;
        setCellDispsPtr_.set(new setCellDisplacements(mesh(), dict));
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
    return rho_;
}


void Foam::solidModel::setCellDisps(fvVectorMatrix& DEqn)
{
    if (setCellDisps().cellIDs().size() > 0)
    {
        DEqn.setValues(setCellDisps().cellIDs(), setCellDisps().cellDisps());
    }
}


void Foam::solidModel::relaxField(volVectorField& D, int iCorr)
{
    if (relaxationMethod_ == "fixed")
    {
        // Fixed under-relaxation
        D.relax();
    }
    else if (relaxationMethod_ == "Aitken")
    {
        // See Aitken method at:
        // http://empire-multiphysics.com/projects/empire/wiki/Aitken_Relaxation
        // and
        // A partitioned solution approach for electro-thermo-
        // problems, Patrick Erbts, Stefan Hartmann, Alexander Duster.

        // Store aitkenResidual previous iteration
        aitkenResidual_.storePrevIter();

        // Calculate new aitkenResidual
        aitkenResidual_ = D.prevIter() - D;

        if (iCorr == 0)
        {
            // Fixed under-relaxation is applied on the first iteration
            aitkenAlpha_ = 1.0;

            if (mesh().relaxField(D.name()))
            {
                aitkenAlpha_ =
                    mesh().fieldRelaxationFactor(D.name());
            }
        }
        else
        {
            const volVectorField aitkenResidualDelta
            (
                aitkenResidual_.prevIter() - aitkenResidual_
            );

            // Update the relaxation factor field
            aitkenAlpha_ =
                aitkenAlpha_*(aitkenResidual_.prevIter() & aitkenResidualDelta)
               /(
                    magSqr(aitkenResidualDelta)
                  + dimensionedScalar("SMALL", dimLength*dimLength, SMALL)
                );

            // Bound alpha between 0.0 and 2.0
            // This may not be necessary but it seems to help convergence
            aitkenAlpha_ = max(0.0, min(2.0, aitkenAlpha_));
        }

        // Relax the field
        D -= aitkenAlpha_*aitkenResidual_;
    }
    else if (relaxationMethod_ == "QuasiNewton")
    {
        // This method is a modified form of the IQNILS by Degroote et al.

        // J. Degroote, K.-J. Bathe and J. Vierendeels.
        // A fluid solid interaction solver with IQN-ILS coupling algorithm.
        // Performance of a new partitioned procedure versus a monolithic
        // procedure in fluid-solid interaction. Computers & Solids

        if (iCorr == 0 || iCorr % QuasiNewtonRestartFreq_ == 0)
        {
            // Clean up data from old time steps

            if (debug)
            {
                Info<< "Modes before clean-up : " << QuasiNewtonT_.size();
            }

            while (true)
            {
                if (QuasiNewtonT_.size())
                {
                    if
                    (
                        runTime().timeIndex() > QuasiNewtonT_[0]
                     || iCorr % QuasiNewtonRestartFreq_ == 0
                    )
                    {
                        for (label i = 0; i < QuasiNewtonT_.size() - 1; i++)
                        {
                            QuasiNewtonT_[i] = QuasiNewtonT_[i + 1];
                            QuasiNewtonV_[i] = QuasiNewtonV_[i + 1];
                            QuasiNewtonW_[i] = QuasiNewtonW_[i + 1];
                        }

                        QuasiNewtonT_.remove();
                        QuasiNewtonV_.remove();
                        QuasiNewtonW_.remove();
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }

            if (debug)
            {
                Info<< ", modes after clean-up : " << QuasiNewtonT_.size()
                    << endl;
            }
        }
        else if (iCorr == 1 || iCorr % QuasiNewtonRestartFreq_ == 1)
        {
            // Set reference in the first coupling iteration
            unrelaxedDRef_ = D;
            DRef_ = D.prevIter();
        }
        else
        {
            // Store the input vector field, defined as the previous iteration
            // D field (after relaxation) minus the Dp previous iteration field
            // in the first iteration (after relaxation)
            QuasiNewtonV_.append
            (
                (D - D.prevIter()) - (unrelaxedDRef_ - DRef_)
            );

            // V should be (from FSI paper):
            // DeltaR^{k-1} = R^{k-1} - R^k = DPrevIt.PrevIt - DPrevIter
            // DeltaR^{k-2} = R^{k-2} - R^k = D.PI.PI.PI - D.PI
            // ...
            // DeltaR^{0} =  R^0 - R^k = DRef - D.prevIter
            // Or in the general paper:
            // V_i = p_k - p_i    for i = 0, 1, ..., k - 1
            // V_{k-1} = p_k - p_{k-1} = D.PI - D.PI.PI
            // V_{k-2} = p_k - p_{k-2} = D.PI - D.PI.PI.PI
            // ...
            // V_{0} = p_k - p_{0} = D.PI - DRef
            // BUT, the implemented code does this:
            // V_i = p_{i+1} - p_0    for i = 0, 1, ..., k - 1
            // V_{k-1} = p_{k} - p_0 = D.PI - DRef
            // V_{k-2} = p_{k-1} - p_0 = D.PI{k-1} - DRef
            // ...
            // V_{0} = p_{1} - p_{0} = D.PI_1 - DRef
            // This means that we just append  the following line each
            // iteration:
            // V_{k-1} = p_{k} - p_0 = D.PI - DRef
            // It this equivalent?
            // We could try implementing it as described in the paper, but this
            // will require D and D.prevIter and their history

            // Store the output vector field, defined as the current iteration
            // D field (before relaxation) minus the D field in the first
            // iteration (before relaxation)
            QuasiNewtonW_.append(D - unrelaxedDRef_);

            // Store the time index
            QuasiNewtonT_.append(runTime().timeIndex());
        }

        if (QuasiNewtonT_.size() > 1)
        {
            // Consider QuasiNewtonV as a matrix V
            // with as columns the items
            // in the DynamicList and calculate the QR-decomposition of V
            // with modified Gram-Schmidt
            label cols = QuasiNewtonV_.size();
            RectangularMatrix<scalar> R(cols, cols, 0.0);
            RectangularMatrix<scalar> C(cols, 1);
            RectangularMatrix<scalar> Rcolsum(1, cols);
            // philipc: do need for dynamic list for Q
            //DynamicList<vectorField> Q(cols);
            List<vectorField> Q(cols);

            for (label i = 0; i < cols; i++)
            {
                //Q.append(QuasiNewtonV_[cols - 1 - i]);
                Q[i] = QuasiNewtonV_[cols - 1 - i];
            }

            for (label i = 0; i < cols; i++)
            {
                // Normalize column i
                R[i][i] = Foam::sqrt(sum(Q[i] & Q[i]));
                Q[i] /= R[i][i];

                // Orthogonalize columns to the right of column i
                for (label j = i+1; j < cols; j++)
                {
                    R[i][j] = sum(Q[i] & Q[j]);
                    Q[j] -= R[i][j]*Q[i];
                }

                // Project minus the residual vector on the Q
                C[i][0] =
                    sum
                    (
                        Q[i]
                      & (
                          D.prevIter().primitiveField()
                        - D.primitiveField()
                        )
                    );
            }

            // Solve the upper triangular system
            for (label j = 0; j < cols; j++)
            {
                Rcolsum[0][j] = 0.0;
                for (label i = 0; i < (j + 1); i++)
                {
                    Rcolsum[0][j] += cmptMag(R[i][j]);
                }
            }
            scalar epsilon = 1.0E-10*max(Rcolsum);
            for (label i = 0; i < cols; i++)
            {
                if (cmptMag(R[i][i]) > epsilon)
                {
                    for (label j = i + 1; j < cols; j++)
                    {
                        R[i][j] /= R[i][i];
                    }
                    C[i][0] /= R[i][i];
                    R[i][i] = 1.0;
                }
            }
            for (label j = (cols - 1); j >= 0; j--)
            {
                if (cmptMag(R[j][j]) > epsilon)
                {
                    for (label i = 0; i < j; i++)
                    {
                        C[i][0] -= C[j][0]*R[i][j];
                    }
                }
                else
                {
                    C[j][0] = 0.0;
                }
            }

            // Update D
            for (label i = 0; i < cols; i++)
            {
                D.primitiveFieldRef() += QuasiNewtonW_[i]*C[cols - 1 - i][0];
            }

            D.correctBoundaryConditions();
        }
        else
        {
            // Fixed under-relaxation during startup
            D.relax();
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::relaxField(volVectorField& D, int iCorr)"
        )   << "relaxationMethod '" << relaxationMethod_ << "' unknown!"
            << " Options are fixed, Aitken or QuasiNewton" << abort(FatalError);
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
        calculatedPointPatchVectorField::typeName
    );
    forAll(D.boundaryField(), patchi)
    {
        if (isA<fixedValueFvPatchVectorField>(D.boundaryField()[patchi]))
        {
            bTypes[patchi] = fixedValuePointPatchVectorField::typeName;
        }
    }
    return bTypes;
}


Foam::dictionary& Foam::solidModel::solidModelDict()
{
    return this->subDict(type_ + "Coeffs");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModel::solidModel
(
    const word& type,
    dynamicFvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonlinear,
    const bool incremental
)
:
//     physicsModel(type, runTime),
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
    thermal_(mesh),
    mechanical_(mesh, nonlinear, incremental),
    Dheader_("D", mesh.time().timeName(), mesh, IOobject::MUST_READ),
    DDheader_("DD", mesh.time().timeName(), mesh, IOobject::MUST_READ),
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
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    DD_
    (
        IOobject
        (
            "DD",
            mesh.time().timeName(),
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
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength/dimTime, vector::zero)
    ),
    pMesh_(mesh),
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
        volPointInterpolation::New(mesh).interpolate(D_)//,
//         pointDBoundaryTypes(D_)
    ),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        volPointInterpolation::New(mesh).interpolate(DD_)//,
//         pointDBoundaryTypes(DD_)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
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
    rho_(thermal_.thermo().rho()),
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
    stabilisationPtr_(),
    solutionTol_
    (
        solidModelDict().lookupOrDefault<scalar>("solutionTolerance", 1e-06)
    ),
    alternativeTol_
    (
        solidModelDict().lookupOrDefault<scalar>("alternativeTolerance", 1e-07)
    ),
    materialTol_
    (
        solidModelDict().lookupOrDefault<scalar>("materialTolerance", 1e-05)
    ),
    infoFrequency_
    (
        solidModelDict().lookupOrDefault<int>("infoFrequency", 100)
    ),
    nCorr_(solidModelDict().lookupOrDefault<int>("nCorrectors", 10000)),
    minCorr_(solidModelDict().lookupOrDefault<int>("minCorrectors", 1)),
    maxIterReached_(0),
    residualFilePtr_(),
    writeResidualField_
    (
        solidModelDict().lookupOrDefault<Switch>("writeResidualField", false)
    ),
    enforceLinear_(false),
    relaxationMethod_
    (
        solidModelDict().lookupOrDefault<word>("relaxationMethod", "fixed")
    ),
    aitkenAlpha_
    (
        IOobject
        (
            "aitkenAlpha",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    ),
    aitkenResidual_
    (
        IOobject
        (
            "aitkenResidual",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    QuasiNewtonRestartFreq_
    (
        solidModelDict().lookupOrDefault<int>("QuasiNewtonRestartFrequency", 25)
    ),
    QuasiNewtonV_(QuasiNewtonRestartFreq_ + 2),
    QuasiNewtonW_(QuasiNewtonRestartFreq_ + 2),
    QuasiNewtonT_(QuasiNewtonRestartFreq_ + 2),
    DRef_
    (
        IOobject
        (
            "DRef",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    unrelaxedDRef_
    (
        IOobject
        (
            "unrelaxedDRef",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    globalPatchesPtrList_()
{
    // Force old time fields to be stored
    D_.oldTime().oldTime();
    DD_.oldTime().oldTime();
    pointD_.oldTime();
    pointDD_.oldTime();
    gradD_.oldTime();
    gradDD_.oldTime();
    sigma_.oldTime();

    // There is an issue where old-old fields are not being written so we will
    // reset the write flag here
    D_.oldTime().oldTime().writeOpt() = IOobject::AUTO_WRITE;
    DD_.oldTime().oldTime().writeOpt() = IOobject::AUTO_WRITE;

    // Print out the relaxation factor
    Info<< "    under-relaxation method: " << relaxationMethod_ << endl;
    if (relaxationMethod_ == "QuasiNewton")
    {
        Info<< "        restart frequency: " << QuasiNewtonRestartFreq_ << endl;
    }

    // If requested, create the residual file
    if (solidModelDict().lookupOrDefault<Switch>("residualFile", false))
    {
        if (Pstream::master())
        {
            Info<< "Creating residual.dat" << endl;
            residualFilePtr_.set
            (
                new OFstream(mesh.time().path()/"residual.dat")
            );
        }
    }

    // Create stabilisation object

    if (!solidModelDict().found("stabilisation"))
    {
        // If the stabilisation sub-dict is not found, we will add it with
        // default settings
        dictionary stabDict;
        stabDict.add("type", "RhieChow");
        stabDict.add("scaleFactor", 0.1);
        solidModelDict().add("stabilisation", stabDict);
    }

    stabilisationPtr_.set
    (
        new momentumStabilisation
        (
            solidModelDict().subDict("stabilisation")
        )
    );

    // If the case is axisymmetric, we will disable solving in the out-of-plane
    // direction
    // PC, 12-Nov-18: disabling the 3rd direction slows down convergence a lot
    // in some elastic cases: disabled for now
    //checkWedges();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidModel::~solidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::solidModel::rho() const
{
    return rho_;
}

const Foam::thermalModel& Foam::solidModel::thermal() const
{
    return thermal_;
}


const Foam::mechanicalModel& Foam::solidModel::mechanical() const
{
    return mechanical_;
}


void Foam::solidModel::DisRequired()
{
    if (!Dheader_.typeHeaderOk<volVectorField>(true))
    {
        FatalErrorIn(type() + "::DisRequired()")
            << "This solidModel requires the 'D' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::DDisRequired()
{
    if (!DDheader_.typeHeaderOk<volVectorField>(true))
    {
        FatalErrorIn(type() + "::DDisRequired()")
            << "This solidModel requires the 'DD' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::makeGlobalPatches
(
    const wordList& patchNames,
    const bool currentConfiguration
) const
{
    globalPatchesPtrList_.setSize(patchNames.size());

    forAll(patchNames, i)
    {
        if (globalPatchesPtrList_.set(i))
        {
            FatalErrorIn
            (
                type() + "::makeGlobalPatches(const wordList&) const"
            )   << "Pointer already set for global patch: "
                << patchNames[i] << "!"
                << abort(FatalError);
        }

        if (currentConfiguration)
        {
            // The global patch will create a standAlone zone based on the
            // current point positions. So we will temporarily move the mesh to
            // the deformed position, then create the globalPatch, then move the
            // mesh back
            const pointField pointsBackup = mesh().points();

            // Lookup patch index
            const label patchID =
                mesh().boundaryMesh().findPatchID(patchNames[i]);
            if (patchID == -1)
            {
                FatalErrorIn("void Foam::solidModel::makeGlobalPatches(...)")
                    << "Patch not found!" << abort(FatalError);
            }

            // Patch point displacement
            const vectorField pointDisplacement
            (
                pointDorPointDD().internalField(),
                mesh().boundaryMesh()[patchID].meshPoints()
            );

            // Calculate deformation point positions
            const pointField newPoints
            (
                mesh().points() + pointDorPointDD().internalField()
            );

            // Move the mesh to deformed position
            // const_cast is justified as it is not our intention to permanently
            // move the mesh; however, it would be better if we did not need it
            mesh().V();
            const_cast<dynamicFvMesh&>(mesh()).movePoints(newPoints);
            const_cast<dynamicFvMesh&>(mesh()).moving(false);
            const_cast<dynamicFvMesh&>(mesh()).setPhi().writeOpt() =
                IOobject::NO_WRITE;

            // Create global patch based on deformed mesh
            globalPatchesPtrList_.set
            (
                i,
                new globalPolyPatch(patchNames[i], mesh())
            );

            // Force creation of standAlonePatch
            globalPatchesPtrList_[i].globalPatch();

            // Move the mesh back
            const_cast<dynamicFvMesh&>(mesh()).movePoints(pointsBackup);
            mesh().V();
            const_cast<dynamicFvMesh&>(mesh()).moving(false);
            const_cast<dynamicFvMesh&>(mesh()).setPhi().writeOpt() =
                IOobject::NO_WRITE;
        }
        else
        {
            globalPatchesPtrList_.set
            (
                i,
                new globalPolyPatch(patchNames[i], mesh())
            );
        }
    }
}


const Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidModel::globalPatches() const
{
    if (globalPatchesPtrList_.empty())
    {
        FatalErrorIn(type() + "::globalPatches() const")
            << "makeGlobalPatches(const wordList&) must be called "
            << "before globalPatch can be called!"
            << abort(FatalError);
    }

    return globalPatchesPtrList_;
}


void Foam::solidModel::clearGlobalPatches() const
{
    globalPatchesPtrList_.clear();
}


Foam::vector Foam::solidModel::pointU(const label pointID) const
{
    pointVectorField pointU
    (
        IOobject
        (
            "pointU",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimVelocity, vector::zero)
    );

    mechanical().volToPoint().interpolate(U(), pointU);

    return pointU.internalField()[pointID];
}


Foam::tmp<Foam::vectorField>
Foam::solidModel::faceZonePointDisplacementIncrement
(
    const label interfaceI
) const
{
    // Create patch point field
    const vectorField patchPointDispIncr
    (
        pointDD().internalField(),
        globalPatches()[interfaceI].patch().meshPoints()
    );

    // Return the global patch field
    return globalPatches()[interfaceI].patchPointToGlobal(patchPointDispIncr);
}


Foam::tmp<Foam::vectorField>
Foam::solidModel::faceZonePointDisplacementOld
(
    const label interfaceI
) const
{
    // Create patch point field
    const vectorField patchPointDispOld
    (
        pointD().oldTime().internalField(),
        globalPatches()[interfaceI].patch().meshPoints()
    );

    // Return the global patch field
    return globalPatches()[interfaceI].patchPointToGlobal(patchPointDispOld);
}


Foam::tmp<Foam::vectorField> Foam::solidModel::faceZoneAcceleration
(
    const label interfaceI
) const
{
    const volVectorField a(fvc::d2dt2(D()));

    return globalPatches()[interfaceI].patchFaceToGlobal
    (
        a.boundaryField()[globalPatches()[interfaceI].patch().index()]
    );
}


void Foam::solidModel::updateTotalFields()
{
    mechanical().updateTotalFields();
}


void Foam::solidModel::end()
{
    if (maxIterReached_ > 0)
    {
        WarningIn(type() + "::end()")
            << "The maximum momentum correctors were reached in "
            << maxIterReached_ << " time-steps" << nl << endl;
    }
    else
    {
        Info<< "The momentum equation converged in all time-steps"
            << nl << endl;
    }
}


Foam::autoPtr<Foam::solidModel> Foam::solidModel::New(dynamicFvMesh& mesh)
{
    word solidModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid model is created, otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary solidProperties
        (
            IOobject
            (
                "solidProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        solidProperties.lookup("solidModel")
            >> solidModelTypeName;
    }

    Info<< "Selecting solidModel " << solidModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solidModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidModel::New(Time&, const word&)"
        )   << "Unknown solidModel type " << solidModelTypeName
            << endl << endl
            << "Valid solidModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidModel>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::solidModel> Foam::solidModel::NewLU(dynamicFvMesh& mesh)
{
    word solidModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid model is created, otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary solidProperties
        (
            IOobject
            (
                "solidProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        solidProperties.lookup("solidModel")
            >> solidModelTypeName;
    }

    Info<< "Selecting solidModel " << solidModelTypeName << endl;

    lagrangianConstructorTable::iterator cstrIter =
        lagrangianConstructorTablePtr_->find(solidModelTypeName);

    if (cstrIter == lagrangianConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidModel::NewLU(Time&, const word&)"
        )   << "Unknown lagrangian solidModel type " << solidModelTypeName
            << endl << endl
            << "Valid lagrangian solidModel types are :" << endl
            << lagrangianConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidModel>(cstrIter()(mesh));
}


void Foam::solidModel::setTraction
(
    fvPatchVectorField& tractionPatch,
    const vectorField& traction
)
{
    if (tractionPatch.type() == solidTractionFvPatchVectorField::typeName)
    {
        solidTractionFvPatchVectorField& patchD =
            refCast<solidTractionFvPatchVectorField>(tractionPatch);

        patchD.traction() = traction;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << tractionPatch.type()
            << " for patch " << tractionPatch.patch().name()
            << " should instead be type "
            << solidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }
}


void Foam::solidModel::setPressure
(
    fvPatchVectorField& pressurePatch,
    const scalarField& pressure
)
{
    if (pressurePatch.type() == solidTractionFvPatchVectorField::typeName)
    {
        solidTractionFvPatchVectorField& patchD =
            refCast<solidTractionFvPatchVectorField>(pressurePatch);

        patchD.pressure() = pressure;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::setPressure\n"
            "(\n"
            "    fvPatchVectorField& pressurePatch,\n"
            "    const vectorField& pressure\n"
            ")"
        )   << "Boundary condition "
            << pressurePatch.type()
            << "for patch" << pressurePatch.patch().name()
            << " should instead be type "
            << solidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }
}


void Foam::solidModel::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    const vectorField patchTraction
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneTraction)
    );

    setTraction(solutionD().boundaryFieldRef()[patchID], patchTraction);
}


void Foam::solidModel::setPressure
(
    const label interfaceI,
    const label patchID,
    const scalarField& faceZonePressure
)
{
    const scalarField patchPressure
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZonePressure)
    );

    setPressure(solutionD().boundaryFieldRef()[patchID], patchPressure);
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const volScalarField& J)
{
    scalar minJ = min(J).value();
    reduce(minJ, minOp<scalar>());

    scalar maxJ = max(J).value();
    reduce(maxJ, maxOp<scalar>());

    if ((minJ < 0.01) || (maxJ > 100))
    {
        Info<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const surfaceScalarField& J)
{
    scalar minJ = min(J).value();
    reduce(minJ, minOp<scalar>());

    scalar maxJ = max(J).value();
    reduce(maxJ, maxOp<scalar>());

    if ((minJ < 0.01) || (maxJ > 100))
    {
        Info<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


void Foam::solidModel::writeFields(const Time& runTime)
{
    // Write strain fields
    // Currently only defined for linear geometry
    if (nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY)
    {
        // Total strain
        volSymmTensorField epsilon("epsilon", symm(gradD()));
        epsilon.write();

        // Equivalent strain
        volScalarField epsilonEq
        (
            "epsilonEq", sqrt((2.0/3.0)*magSqr(dev(epsilon)))
        );
        epsilonEq.write();

        Info<< "Max epsilonEq = " << gMax(epsilonEq) << endl;
    }

    // Calculate equivalent (von Mises) stress
    volScalarField sigmaEq
    (
        "sigmaEq", sqrt((3.0/2.0)*magSqr(dev(sigma())))
    );
    sigmaEq.write();

    Info<< "Max sigmaEq (von Mises stress) = " << gMax(sigmaEq) << endl;

    // If asked, write the residual field
    if (writeResidualField_)
    {
        const volVectorField& D = solutionD();
        scalar denom =
            gMax(mag(D.primitiveField() - D.oldTime().primitiveField()));
        if (denom < SMALL)
        {
            denom = max(gMax(mag(D.primitiveField())), SMALL);
        }

        const volVectorField residualD
        (
            "residualD",
            (D - D.prevIter())/denom
        );

        Info<< "Writing residualD field" << endl;
        residualD.write();
    }

//     physicsModel::writeFields(runTime);
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

            if
            (
                returnReduce(mesh().boundaryMesh()[patchI].size(), sumOp<int>())
            )
            {
                continue;
            }

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());
            const vector k(0, 0, 1);

            if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].z() = 0;
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

    twoDPointCorrector twoDCorrector(mesh());
    twoDCorrector.correctPoints(newPoints);
    twoDCorrector.correctPoints(pointDDI);
    mesh().movePoints(newPoints);
    mesh().V00();
    mesh().moving(false);
    mesh().setPhi().writeOpt() = IOobject::NO_WRITE;
}


const Foam::dictionary& Foam::solidModel::solidModelDict() const
{
    return this->subDict(type_ + "Coeffs");
}

// ************************************************************************* //
