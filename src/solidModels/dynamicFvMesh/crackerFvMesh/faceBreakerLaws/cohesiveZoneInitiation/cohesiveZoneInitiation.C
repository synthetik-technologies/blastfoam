/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cohesiveZoneInitiation.H"
#include "addToRunTimeSelectionTable.H"
#include "crackerFvMesh.H"
#include "solidCohesiveFvPatchVectorField.H"
#include "cohesivePolyPatch.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cohesiveZoneInitiation, 0);
    addToRunTimeSelectionTable
    (
        faceBreakerLaw, cohesiveZoneInitiation, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::cohesiveZoneInitiation::calcCohesivePatchID() const
{
    if (cohesivePatchIDPtr_)
    {
        FatalErrorIn
        (
            "void Foam::cohesiveZoneInitiation::calcCohesivePatchID() const"
        ) << "pointer already set" << abort(FatalError);
    }

    const fvMesh& mesh = this->mesh();

    cohesivePatchIDPtr_ = new label(-1);
    label& cohesivePatchID = *cohesivePatchIDPtr_;

    forAll (mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].type() == cohesivePolyPatch::typeName)
        {
            cohesivePatchID = patchI;
            break;
        }
    }

    if (cohesivePatchID == -1)
    {
        FatalErrorIn
        (
            "void Foam::cohesiveZoneInitiation::calcCohesivePatchID() const"
        )   << "boundary patch of type cohesive not found" << abort(FatalError);
    }
}


Foam::label Foam::cohesiveZoneInitiation::cohesivePatchID() const
{
    if (!cohesivePatchIDPtr_)
    {
        calcCohesivePatchID();
    }

    return *cohesivePatchIDPtr_;
}


const Foam::cohesiveZoneModel&
Foam::cohesiveZoneInitiation::cohesiveZone() const
{
    // Const reference to the mesh
    const fvMesh& mesh = this->mesh();

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(mesh);

    // Lookup displacement field
    if (solMod.incremental())
    {
        return
            refCast<const solidCohesiveFvPatchVectorField>
            (
                mesh.lookupObject<volVectorField>
                (
                    "DD"
                ).boundaryField()[cohesivePatchID()]
            ).cohesiveZone();
    }
    else
    {
        return
            refCast<const solidCohesiveFvPatchVectorField>
            (
                mesh.lookupObject<volVectorField>
                (
                    "D"
                ).boundaryField()[cohesivePatchID()]
            ).cohesiveZone();
    }
}


void Foam::cohesiveZoneInitiation::calcAllFacesToBreak() const
{
    if (facesToBreakPtr_ || coupledFacesToBreakPtr_)
    {
        FatalErrorIn
        (
            "void Foam::cohesiveZoneInitiation::calcAllFacesToBreak() const"
        ) << "pointer already set" << abort(FatalError);
    }

    // First, we check if any internal faces need to be broken, then we will
    // check if any coupled (processor boundary) faces need to be broken.
    // If we find more than one face with a traction fraction greater than 1.0,
    // we will select the face with the highest value.
    // Note: the traction fraction is the traction divided by the strength of
    // the face, so a value greater than or equal to 1.0 indicates that the face
    // should break.

    int nFacesToBreak = 0;
    int nCoupledFacesToBreak = 0;

    // Cast the mesh to a crackerFvMesh

    if (!isA<crackerFvMesh>(this->mesh()))
    {
        FatalErrorIn("Foam::label Foam::cohesiveZoneInitiation::updateMesh()")
            << "Mesh should be of type: " << crackerFvMesh::typeName
            << abort(FatalError);
    }

    const crackerFvMesh& mesh =
        dynamicCast<const crackerFvMesh>(this->mesh());

    // Face unit normals
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Get traction fraction from cohesive zone model
    cohesiveZone().updateMeshTraction();
    surfaceScalarField tracFrac
    (
            cohesiveZone().initiationTractionFraction()
    );
    const surfaceVectorField& traction = cohesiveZone().meshTraction();

    // Apply crack path limiting if specified
    if (pathLimiterPtr_.valid())
    {
        pathLimiterPtr_().clearOut();
        tracFrac *= pathLimiterPtr_().facesAllowedToBreak();
    }

    const scalar maxTracFrac = gMax(tracFrac.internalField());

    Info<< nl << "Max traction fraction: " << maxTracFrac << endl;

    label faceToBreakIndex = -1;
    scalar faceToBreakTracFrac = 0.0;

    if (maxTracFrac > 1.0)
    {
        // Find face with maximum traction fraction

        forAll(tracFrac, faceI)
        {
            if (mag(tracFrac[faceI] - maxTracFrac) < SMALL)
            {
                faceToBreakIndex = faceI;
                faceToBreakTracFrac = maxTracFrac;
                break;
            }
        }

        if (faceToBreakIndex != -1)
        {
            Pout<< "    faceToBreakIndex: " << faceToBreakIndex << nl
                << "    faceToBreakLocation: " << mesh.Cf()[faceToBreakIndex]
                << nl
                << "    faceToBreakEffTracFrac: " << maxTracFrac
                << endl;
        }

        nFacesToBreak = 1;

        if (Pstream::parRun())
        {
            // Find processor with greatest traction fraction

            bool procHasFaceToBreak = false;

            if (mag(faceToBreakTracFrac - maxTracFrac) < SMALL)
            {
                procHasFaceToBreak = true;
            }

            // Check if maximum is present on more then one processors

            label procID = Pstream::nProcs();

            if (procHasFaceToBreak)
            {
                procID = Pstream::myProcNo();
            }

            const label minProcID = returnReduce<label>(procID, minOp<label>());

            if (procID != minProcID)
            {
                nFacesToBreak = 0;
            }
        }
    }

    // Check coupled (processor) patches

    label coupledFaceToBreakIndex = -1;
    scalar coupledFaceToBreakTracFrac = 0.0;
    scalar maxCoupledTracFrac = 0.0;

    if (allowCoupledFaces_ && Pstream::parRun())
    {
        forAll(mesh.boundary(), patchI)
        {
            // Find coupled face with the largest traction fraction
            if (mesh.boundary()[patchI].coupled())
            {
                const scalarField& pTracFrac =
                    tracFrac.boundaryField()[patchI];
                const label start = mesh.boundaryMesh()[patchI].start();

                forAll(pTracFrac, faceI)
                {
                    if
                    (
                        pTracFrac[faceI]
                      > coupledFaceToBreakTracFrac
                    )
                    {
                        coupledFaceToBreakIndex = faceI + start;
                        coupledFaceToBreakTracFrac =
                            pTracFrac[faceI];
                    }
                }
            }
        }

        // Get global coupled traction fraction
        maxCoupledTracFrac =
            returnReduce(coupledFaceToBreakTracFrac, maxOp<scalar>());

        Info<< "Max coupled traction fraction: " << maxCoupledTracFrac
            << endl;

        if (maxCoupledTracFrac > 1.0)
        {
            nCoupledFacesToBreak = 1;

            // Find processor with the greatest traction fraction on a processor
            // face

            bool procHasCoupledFaceToBreak = false;

            if
            (
                mag(maxCoupledTracFrac - coupledFaceToBreakTracFrac)
              < SMALL
            )
            {
                procHasCoupledFaceToBreak = true;
            }

            // Check if maximum is present on more then one processors

            label procID = Pstream::nProcs();

            if (procHasCoupledFaceToBreak)
            {
                procID = Pstream::myProcNo();
            }

            label minProcID = returnReduce<label>(procID, minOp<label>());

            if (procID != minProcID)
            {
                nCoupledFacesToBreak = 0;
            }
            // Now procID contains the minimum ID of the processor with the
            // coupled face to break on it
        }
    }

    // If both an internal and coupled face want to break then we pick the face
    // with the greatest traction fraction

    if (maxCoupledTracFrac > maxTracFrac)
    {
        // Break coupled face
        nFacesToBreak = 0;
    }
    else
    {
        // Break internal face
        nCoupledFacesToBreak = 0;
    }

    // Make sure that coupled faces are broken in pairs

    // ngbProc will store the neighbor processor number for the current
    // processor if the coupled face to break is on the current processor
    // All the other indices of the ngbProc list will remain -1
    labelList ngbProc(Pstream::nProcs(), -1);
    // index will store the local faceID of the coupledFace to break on the
    // current processor
    labelList index(Pstream::nProcs(), -1);

    // This if statement will only be entered if we want to break the coupled
    // face and if the coupled face is on the current processor
    if (nCoupledFacesToBreak)
    {
        label patchID =
            mesh.boundaryMesh().whichPatch(coupledFaceToBreakIndex);

        if (patchID == -1)
        {
            FatalErrorIn("cohesiveZoneInitiation::calcAllFacesToBreak()")
                << "something is wrong: patchID is -1 for coupled face"
                << abort(FatalError);
        }

        label start = mesh.boundaryMesh()[patchID].start();
        label localIndex = coupledFaceToBreakIndex - start;

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchID]);

        label ngbProcNo = procPatch.neighbProcNo();

        ngbProc[Pstream::myProcNo()] = ngbProcNo;
        index[Pstream::myProcNo()] = localIndex;
    }

    // Enter this if condition for all processor if we are planning to break a
    // coupled face in this iteration
    if (returnReduce(nCoupledFacesToBreak, maxOp<label>()))
    {
        // syncing the ngbProc and index lists so that all the processors have
        // them
        reduce(ngbProc, maxOp<labelList>());
        reduce(index, maxOp<labelList>());

        // Going through all the processors to find the processor that has the
        // coupled part of the coupled face to break in this iteration
        // This is done by comparing the current proc no with the processor no
        // stored in the ngbProc list and then finding the corresponding coupled
        // face on the neighbour processor
        // The assumption here is that the faces on both processor boundaries
        // have the same local faceID
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                if (ngbProc[procI] == Pstream::myProcNo())
                {
                    forAll(mesh.boundaryMesh(), patchI)
                    {
                        if
                        (
                            mesh.boundaryMesh()[patchI].type()
                         == processorPolyPatch::typeName
                        )
                        {
                            const processorPolyPatch& procPatch =
                                refCast<const processorPolyPatch>
                                (
                                    mesh.boundaryMesh()[patchI]
                                );

                            label ngbProcNo = procPatch.neighbProcNo();

                            if (ngbProcNo == procI)
                            {
                                label start =
                                    mesh.boundaryMesh()[patchI].start();
                                coupledFaceToBreakIndex = start + index[procI];
                                nCoupledFacesToBreak = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    // Note: It is not necessary to scale the tractions as the cohesive boundary
    // condition will do this

    if (returnReduce(nCoupledFacesToBreak, sumOp<label>()) > 2)
    {
        FatalErrorIn("cohesiveZoneInitiation::calcAllFacesToBreak()")
            << "More than two processors are trying to break a coupled face"
            << abort(FatalError);
    }

    // Now we return the faces to break

    if (nFacesToBreak)
    {
        facesToBreakPtr_ = new labelList(nFacesToBreak, faceToBreakIndex);

        coupledFacesToBreakPtr_ = new labelList(0);

        // Store tractions of face to break

        facesToBreakTractionsPtr_ =
            new List<vector>
            (
                nFacesToBreak, traction.internalField()[faceToBreakIndex]
            );

        coupledFacesToBreakTractionsPtr_ = new List<vector>(0);

        // Store normals of the face to break

        facesToBreakNormalsPtr_ =
            new List<vector>
            (
                nFacesToBreak, n.internalField()[faceToBreakIndex]
            );

        coupledFacesToBreakNormalsPtr_ = new List<vector>(0);
    }
    else if (nCoupledFacesToBreak)
    {
        facesToBreakPtr_ = new labelList(0);

        coupledFacesToBreakPtr_ =
            new labelList(nCoupledFacesToBreak, coupledFaceToBreakIndex);

        // Store tractions of the face to break

        facesToBreakTractionsPtr_ = new List<vector>(0);

        const label patchID =
            mesh.boundaryMesh().whichPatch(coupledFaceToBreakIndex);
        const label start = mesh.boundaryMesh()[patchID].start();

        coupledFacesToBreakTractionsPtr_ =
            new List<vector>
            (
                nCoupledFacesToBreak,
                traction.boundaryField()[patchID]
                [
                    coupledFaceToBreakIndex - start
                ]
            );

        // Store normals of the face to break

        facesToBreakNormalsPtr_ = new List<vector>(0);

        coupledFacesToBreakNormalsPtr_ =
            new List<vector>
            (
                nCoupledFacesToBreak,
                n.boundaryField()[patchID][coupledFaceToBreakIndex - start]
            );
    }
    else
    {
        facesToBreakPtr_ = new labelList(0);
        coupledFacesToBreakPtr_ = new labelList(0);
        facesToBreakTractionsPtr_ = new List<vector>(0);
        coupledFacesToBreakTractionsPtr_ = new List<vector>(0);
        facesToBreakNormalsPtr_ = new List<vector>(0);
        coupledFacesToBreakNormalsPtr_ = new List<vector>(0);
    }

    if (debug)
    {
        const labelList& f = *facesToBreakPtr_;
        forAll(f, faceI)
        {
            Pout<< "coupled face to break: " << mesh.Cf()[f[faceI]] << endl;
        }

        const labelList& cf = *coupledFacesToBreakPtr_;

        forAll(cf, faceI)
        {
            const label patchID = mesh.boundaryMesh().whichPatch(cf[faceI]);
            const label start = mesh.boundaryMesh()[patchID].start();

            Pout<< "coupled face to break: "
                << mesh.boundaryMesh()[patchID].faceCentres()[cf[faceI] - start]
                << endl;
        }
    }
}


void Foam::cohesiveZoneInitiation::clearOut()
{
    deleteDemandDrivenData(cohesivePatchIDPtr_);
    deleteDemandDrivenData(facesToBreakPtr_);
    deleteDemandDrivenData(coupledFacesToBreakPtr_);
    deleteDemandDrivenData(facesToBreakTractionsPtr_);
    deleteDemandDrivenData(coupledFacesToBreakTractionsPtr_);
    deleteDemandDrivenData(facesToBreakNormalsPtr_);
    deleteDemandDrivenData(coupledFacesToBreakNormalsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from dictionary
Foam::cohesiveZoneInitiation::cohesiveZoneInitiation
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    faceBreakerLaw(name, mesh, dict),
    cohesivePatchIDPtr_(NULL),
    allowCoupledFaces_(dict.lookupOrDefault<Switch>("allowCoupledFaces", true)),
    facesToBreakPtr_(NULL),
    coupledFacesToBreakPtr_(NULL),
    facesToBreakTractionsPtr_(NULL),
    coupledFacesToBreakTractionsPtr_(NULL),
    facesToBreakNormalsPtr_(NULL),
    coupledFacesToBreakNormalsPtr_(NULL),
    pathLimiterPtr_(NULL)
{
    if (!allowCoupledFaces_)
    {
        WarningIn("cohesiveZoneInitiation::cohesiveZoneInitiation(...)")
            << name << ": allowCoupledFaces is false" << endl;
    }

    // If specified, create crackPathLimiter
    if (dict.found("crackPathLimiter"))
    {
        pathLimiterPtr_ =
            crackPathLimiter::New
            (
                "law", mesh, dict.subDict("crackPathLimiter")
            );
    }
    else
    {
        Info<< "crackPathLimiter not specified" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::cohesiveZoneInitiation::~cohesiveZoneInitiation()
{
    clearOut();
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //


const Foam::labelList& Foam::cohesiveZoneInitiation::facesToBreak() const
{
    if (!facesToBreakPtr_)
    {
        calcAllFacesToBreak();
    }

    return *facesToBreakPtr_;
}


const Foam::labelList& Foam::cohesiveZoneInitiation::coupledFacesToBreak() const
{
    if (!coupledFacesToBreakPtr_)
    {
        calcAllFacesToBreak();
    }

    return *coupledFacesToBreakPtr_;
}


const Foam::List<Foam::vector>&
Foam::cohesiveZoneInitiation::facesToBreakTractions() const
{
    if (!facesToBreakTractionsPtr_)
    {
        calcAllFacesToBreak();
    }

    return *facesToBreakTractionsPtr_;
}


const Foam::List<Foam::vector>&
Foam::cohesiveZoneInitiation::coupledFacesToBreakTractions() const
{
    if (!coupledFacesToBreakTractionsPtr_)
    {
        calcAllFacesToBreak();
    }

    return *coupledFacesToBreakTractionsPtr_;
}


const Foam::List<Foam::vector>&
Foam::cohesiveZoneInitiation::facesToBreakNormals() const
{
    if (!facesToBreakNormalsPtr_)
    {
        calcAllFacesToBreak();
    }

    return *facesToBreakNormalsPtr_;
}


const Foam::List<Foam::vector>&
Foam::cohesiveZoneInitiation::coupledFacesToBreakNormals() const
{
    if (!coupledFacesToBreakNormalsPtr_)
    {
        calcAllFacesToBreak();
    }

    return *coupledFacesToBreakNormalsPtr_;
}


// ************************************************************************* //
