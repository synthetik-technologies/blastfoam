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
#include "wedgeFvPatchFields.H"
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::solidModel::checkWedges(const volVectorField& D)
{
    const fvMesh& mesh = D.mesh();
    label nWedgePatches = 0;
    vector wedgeDirVec = vector::zero;

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        if (isA<wedgePolyPatch>(pp))
        {
            const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>(pp);

            nWedgePatches++;
            wedgeDirVec += cmptMag(wpp.centreNormal());

            // Make sure that solidWedge is used instead of wedge
            if (isA<wedgeFvPatchVectorField>(D.boundaryField()[patchI]))
            {
                WarningInFunction
                    << "solidWedge should be used on displacement solution "
                    << "field wedge patches as non-orthogonal corrections "
                    << "are important!" << endl;
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
                FatalErrorInFunction
                    << "There are no faces on the wedge patch "
                    << mesh.boundaryMesh()[patchI].name()
                    << " on this processor:"
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


Foam::wordList Foam::solidModel::pointDBoundaryTypes
(
    const volVectorField& D
)
{
    wordList bTypes
    (
        D.boundaryField().size(),
        "calculated"
    );
    forAll(D.boundaryField(), patchi)
    {
        if (D.boundaryField()[patchi].fixesValue())
        {
            bTypes[patchi] = "fixedValue";
        }
    }
    return bTypes;
}


// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

Foam::thermalModel& Foam::solidModel::thermal()
{
    return thermal_;
}


Foam::volScalarField& Foam::solidModel::rho()
{
    return thermal_.rho();
}


Foam::dictionary& Foam::solidModel::solidModelDict()
{
    if (this->isDict(type_ + "Coeffs"))
    {
        return this->subDict(type_ + "Coeffs");
    }
    return *this;
}

void Foam::solidModel::readDict()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModel::solidModel
(
    const word& type,
    fvMesh& mesh,
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
            IOobject::AUTO_WRITE
        )
    ),
    mesh_(mesh),
    type_(type),
    thermal_(mesh, isSolid),
    globalPatches_(globalPolyBoundaryMesh::New(mesh))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidModel::~solidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidModel::isRequired
(
    const volVectorField& D,
    const word& type
) const
{
    if (!const_cast<volVectorField&>(D).headerOk())
    {
        FatalErrorInFunction
            << type << " requires the 'D' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::initialize()
{}

const Foam::volScalarField& Foam::solidModel::rho() const
{
    return thermal_.rho();
}

const Foam::thermalModel& Foam::solidModel::thermal() const
{
    return thermal_;
}

const Foam::dictionary& Foam::solidModel::solidModelDict() const
{
    return this->optionalSubDict(type_ + "Coeffs");
}


void Foam::solidModel::updateTotalFields()
{
    thermal().correct();
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


Foam::scalar Foam::solidModel::newDeltaT() const
{
    return great;
}

void Foam::solidModel::moveMesh
(
    const pointField& oldPoints,
    const volVectorField& DD,
    pointVectorField& pointDD
)
{}


bool Foam::solidModel::checkEnforceLinear(const volScalarField& J)
{
    scalar minJ = min(J).value();
    scalar maxJ = max(J).value();
    if ((minJ < 0.01) || (maxJ > 100))
    {
        DebugInfo<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear(true);
    }

    return enforceLinear();
}


bool Foam::solidModel::checkEnforceLinear(const surfaceScalarField& J)
{
    scalar minJ = min(J).value();
    scalar maxJ = max(J).value();
    if ((minJ < 0.01) || (maxJ > 100))
    {
        DebugInfo<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear(true);
    }

    return enforceLinear();
}


// ************************************************************************* //
