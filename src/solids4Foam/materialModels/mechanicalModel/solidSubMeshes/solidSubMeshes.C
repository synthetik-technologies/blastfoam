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

#include "solidSubMeshes.H"
#include "twoDPointCorrector.H"
#include "wedgePolyPatch.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::solidSubMeshes::makeSubMeshes() const
{
    if (!subsetMeshes_.empty())
    {
        FatalErrorInFunction
            << "sub-meshes already exist" << abort(FatalError);
    }

    if (cellZoneNames_.size() == 1)
    {
        FatalErrorInFunction
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    // Check that each cell is in exactly one cellZone
    checkCellZones();

    labelList region(baseMesh().nCells(), -1);

    const meshCellZones& cellZones = baseMesh().cellZones();
    forAll(cellZoneNames_, matI)
    {
        const label cellZoneID = cellZones.findIndex(cellZoneNames_[matI]);

        if (cellZoneID < 0)
        {
            FatalErrorInFunction
                << "cellZone not found for material " << cellZoneNames_[matI]
                << abort(FatalError);
        }

        const labelList& curCellZone =
            baseMesh().cellZones()[cellZoneID];

        forAll(curCellZone, cI)
        {
            region[curCellZone[cI]] = matI;
        }
    }

    subsetMeshes_.setSize(cellZoneNames_.size());

    forAll(subsetMeshes_, matI)
    {
        word subMeshName = cellZoneNames_[matI];

        if (Pstream::parRun())
        {
            subMeshName =
                "proc" + Foam::name(Pstream::myProcNo()) + "_"
               + cellZoneNames_[matI];
        }

        subsetMeshes_.set
        (
            matI,
            new meshSubset(baseMesh())
        );

        subsetMeshes_[matI].setLargeCellSubset(region, matI);
        dynamicCast<regIOobject>
        (
            dynamicCast<polyMesh>(subsetMeshes_[matI].subMesh())
        ).rename(cellZoneNames_[matI]);

        {
            // Try accessing dictionaries
            objectRegistry obr
            (
                IOobject
                (
                    cellZoneNames_[matI],
                    baseMesh().time().timeName(),
                    baseMesh().time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            );

            IOdictionary fvSchemesDict
            (
                IOobject
                (
                    "fvSchemes",
                    obr.time().system(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            );

            if (!fvSchemesDict.headerOk())
            {
                Info << "Cannot read " << fvSchemesDict.path()
                    << ".  Copy from base" << endl;

                IOdictionary fvSchemesBase
                (
                    IOobject
                    (
                        "fvSchemes",
                        baseMesh().time().system(),
                        baseMesh().time(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                fvSchemesDict = fvSchemesBase;
                fvSchemesDict.regIOobject::write();
            }
            dynamicCast<fvSchemes>
            (
                const_cast<fvMesh&>(subsetMeshes_[matI].subMesh())
            ).readOpt() = IOobject::READ_IF_PRESENT;
            dynamicCast<fvSchemes>
            (
                const_cast<fvMesh&>(subsetMeshes_[matI].subMesh())
            ).read();

            IOdictionary fvSolutionDict
            (
                IOobject
                (
                    "fvSolution",
                    obr.time().system(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            );

            if (!fvSolutionDict.headerOk())
            {
                Info << "Cannot read " << fvSolutionDict.path()
                    << ".  Copy from base" << endl;

                IOdictionary fvSolutionBase
                (
                    IOobject
                    (
                        "fvSolution",
                        baseMesh().time().system(),
                        baseMesh().time(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                fvSolutionDict = fvSolutionBase;
                fvSolutionDict.regIOobject::write();
            }
            dynamicCast<fvSolution>
            (
                const_cast<fvMesh&>(subsetMeshes_[matI].subMesh())
            ).readOpt() = IOobject::READ_IF_PRESENT;
            dynamicCast<fvSolution>
            (
                const_cast<fvMesh&>(subsetMeshes_[matI].subMesh())
            ).read();
        }
    }
}


void Foam::solidSubMeshes::makeSubMeshVolToPoint() const
{
    if (!subMeshVolToPoint_.empty())
    {
        FatalErrorInFunction
            << "sub-meshes already exist" << abort(FatalError);
    }

    if (cellZoneNames_.size() == 1)
    {
        FatalErrorInFunction
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshVolToPoint_.setSize(cellZoneNames_.size());

    forAll(subMeshVolToPoint_, matI)
    {
        subMeshVolToPoint_.set
        (
            matI,
            new volPointInterp(subMeshes()[matI].subMesh())
        );
    }
}


void Foam::solidSubMeshes::checkCellZones() const
{
    if (cellZoneNames_.size() == 1)
    {
        // Cell zones need not be defined if there is only one material
        return;
    }

    // We will check that every cell is in exactly one cellZone
    volScalarField nCellZones
    (
        IOobject
        (
            "nCellZones",
            baseMesh().time().timeName(),
            baseMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        baseMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    );

    scalarField& nCellZonesI = nCellZones.primitiveFieldRef();

    const meshCellZones& cellZones = baseMesh().cellZones();
    forAll(cellZoneNames_, matI)
    {
        const label cellZoneID =
            cellZones.findIndex(cellZoneNames_[matI]);

        if (cellZoneID < 0)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneNames_[matI]
                << " not found for material " << cellZoneNames_[matI]
                << abort(FatalError);
        }

        const labelList& curCellZone =
            baseMesh().cellZones()[cellZoneID];

        forAll(curCellZone, cI)
        {
            const label cellID = curCellZone[cI];
            nCellZonesI[cellID] = nCellZonesI[cellID] + 1.0;
        }
    }

    if (mag(gMin(nCellZonesI)) < small)
    {
        FatalErrorInFunction
            << "There are cells that are not in a material cellZone!"
            << abort(FatalError);
    }

    if (mag(gMax(nCellZonesI) - 1) > small)
    {
        FatalErrorInFunction
            << "There are cells that are in more than one material cellZone!"
            << abort(FatalError);
    }
}


void Foam::solidSubMeshes::calcSubMeshSigma() const
{
    if (!subMeshSigma_.empty())
    {
        FatalErrorInFunction
            << "pointer list already set" << abort(FatalError);
    }

    if (cellZoneNames_.size() == 1)
    {
        // Sub-meshes should not be defined if there is only one material
        FatalErrorInFunction
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshSigma_.setSize(cellZoneNames_.size());

    const meshSubsetList& subMeshes = this->subMeshes();

    forAll(cellZoneNames_, matI)
    {
        subMeshSigma_.set
        (
            matI,
            new volSymmTensorField
            (
                IOobject
                (
                    "sigma",
                    subMeshes[matI].subMesh().time().timeName(),
                    subMeshes[matI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[matI].subMesh(),
                dimensionedSymmTensor
                (
                    "zero",
                    dimForce/dimArea,
                    symmTensor::zero
                )
            )
        );
    }
}


void Foam::solidSubMeshes::calcSubMeshSigmaf() const
{
    if (!subMeshSigmaf_.empty())
    {
        FatalErrorInFunction
            << "pointer list already set" << abort(FatalError);
    }

    if (cellZoneNames_.size() == 1)
    {
        // Sub-meshes should not be defined if there is only one material
        FatalErrorInFunction
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshSigmaf_.setSize(cellZoneNames_.size());

    const meshSubsetList& subMeshes = this->subMeshes();

    forAll(cellZoneNames_, matI)
    {
        subMeshSigmaf_.set
        (
            matI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "sigmaf",
                    subMeshes[matI].subMesh().time().timeName(),
                    subMeshes[matI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[matI].subMesh(),
                dimensionedSymmTensor
                (
                    "zero",
                    dimForce/dimArea,
                    symmTensor::zero
                )
            )
        );
    }
}


void Foam::solidSubMeshes::calcSubMeshD() const
{
    if (!subMeshD_.empty())
    {
        FatalErrorInFunction
            << "pointer list already set" << abort(FatalError);
    }

    if (cellZoneNames_.size() == 1)
    {
        FatalErrorInFunction
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshD_.setSize(cellZoneNames_.size());

    const meshSubsetList& subMeshes = this->subMeshes();

    // The subMeshD field can represent D or DD
    word Dname = "D";
    if (incremental_)
    {
        Dname = "DD";
    }

    forAll(cellZoneNames_, matI)
    {
        subMeshD_.set
        (
            matI,
            new volVectorField
            (
                IOobject
                (
                    Dname,
                    subMeshes[matI].subMesh().time().timeName(),
                    subMeshes[matI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[matI].subMesh(),
                dimensionedVector("zero", dimLength, vector::zero)
            )
        );
    }
}


void Foam::solidSubMeshes::calcSubMeshGradD() const
{
    if (!subMeshGradD_.empty())
    {
        FatalErrorInFunction
            << "pointer list already set" << abort(FatalError);
    }

    if (cellZoneNames_.size() == 1)
    {
        FatalErrorInFunction
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshGradD_.setSize(cellZoneNames_.size());

    const meshSubsetList& subMeshes = this->subMeshes();

    // The subMeshD field can represent D or DD
    word gradDname = "grad(D)";
    if (incremental_)
    {
        gradDname = "grad(DD)";
    }

    forAll(cellZoneNames_, matI)
    {
        subMeshGradD_.set
        (
            matI,
            new volTensorField
            (
                IOobject
                (
                    gradDname,
                    subMeshes[matI].subMesh().time().timeName(),
                    subMeshes[matI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[matI].subMesh(),
                dimensionedTensor("zero", dimless, tensor::zero)
            )
        );
    }
}


void Foam::solidSubMeshes::calcSubMeshGradDf() const
{
    if (!subMeshGradDf_.empty())
    {
        FatalErrorInFunction
            << "pointer list already set" << abort(FatalError);
    }

    if (cellZoneNames_.size() == 1)
    {
        FatalErrorInFunction
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshGradDf_.setSize(cellZoneNames_.size());

    const meshSubsetList& subMeshes = this->subMeshes();

    // The subMeshD field can represent D or DD
    word gradDname = "grad(D)f";
    if (incremental_)
    {
        gradDname = "grad(DD)f";
    }

    forAll(cellZoneNames_, matI)
    {
        subMeshGradDf_.set
        (
            matI,
            new surfaceTensorField
            (
                IOobject
                (
                    gradDname,
                    subMeshes[matI].subMesh().time().timeName(),
                    subMeshes[matI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes[matI].subMesh(),
                dimensionedTensor("zero", dimless, tensor::zero)
            )
        );
    }
}


void Foam::solidSubMeshes::calcSubMeshPointD() const
{
    if (!subMeshPointD_.empty())
    {
        FatalErrorInFunction
            << "pointer list already set" << abort(FatalError);
    }

    if (cellZoneNames_.size() == 1)
    {
        FatalErrorInFunction
            << "There should be no need for subMeshes when there is only one "
            << "material" << abort(FatalError);
    }

    subMeshPointD_.setSize(cellZoneNames_.size());

    const meshSubsetList& subMeshes = this->subMeshes();

    // The subMeshD field can represent D or DD
    word pointDname = "pointD";
    if (incremental_)
    {
        pointDname = "pointDD";
    }

    forAll(cellZoneNames_, matI)
    {
        subMeshPointD_.set
        (
            matI,
            new pointVectorField
            (
                IOobject
                (
                    pointDname,
                    subMeshes[matI].subMesh().time().timeName(),
                    subMeshes[matI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pointMesh::New(subMeshes[matI].subMesh()),
                dimensionedVector("zero", dimLength, vector::zero)
            )
        );
    }
}


void Foam::solidSubMeshes::calcInterfaceShadowIDs() const
{
    if
    (
        !interfaceShadowSubMeshID_.empty()
     || !interfaceShadowPatchID_.empty()
     || !interfaceShadowFaceID_.empty()
    )
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    interfaceShadowSubMeshID_.setSize(cellZoneNames_.size());
    interfaceShadowPatchID_.setSize(cellZoneNames_.size());
    interfaceShadowFaceID_.setSize(cellZoneNames_.size());

    // Reverse maps for the interface faces
    labelList baseMeshShadSubMeshIDs = labelList(baseMesh().nInternalFaces(), -1);
    labelList baseMeshShadPatchIDs = labelList(baseMesh().nInternalFaces(), -1);
    labelList baseMeshShadFaceIDs = labelList(baseMesh().nInternalFaces(), -1);

    forAll(cellZoneNames_, matI)
    {
        const fvMesh& subMesh = subMeshes()[matI].subMesh();

        const labelList& patchMap = subMeshes()[matI].patchMap();
        const labelList& faceMap = subMeshes()[matI].faceMap();

        bool uniquePatchFound = false;

        forAll(subMesh.boundaryMesh(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                if (uniquePatchFound)
                {
                    FatalErrorInFunction
                        << "There are more than one interface patches!"
                        << abort(FatalError);
                }

                uniquePatchFound = true;

                interfaceShadowSubMeshID_.set
                (
                    matI,
                    new labelList(subMesh.boundaryMesh()[patchI].size(), -1)
                );

                interfaceShadowPatchID_.set
                (
                    matI,
                    new labelList(subMesh.boundaryMesh()[patchI].size(), -1)
                );

                interfaceShadowFaceID_.set
                (
                    matI,
                    new labelList(subMesh.boundaryMesh()[patchI].size(), -1)
                );

                labelList& shadSubMeshIDs = interfaceShadowSubMeshID_[matI];
                labelList& shadPatchIDs = interfaceShadowPatchID_[matI];
                labelList& shadFaceIDs = interfaceShadowFaceID_[matI];

                const label start = subMesh.boundaryMesh()[patchI].start();

                forAll(shadSubMeshIDs, faceI)
                {
                    const label baseFaceID = faceMap[start + faceI];

                    if (baseMesh().isInternalFace(baseFaceID))
                    {
                        // Check if the face has been set in the baseMesh lists
                        if (baseMeshShadSubMeshIDs[baseFaceID] == -1)
                        {
                            // Store the local IDs in the baseMesh lists
                            baseMeshShadSubMeshIDs[baseFaceID] = matI;
                            baseMeshShadPatchIDs[baseFaceID] = patchI;
                            baseMeshShadFaceIDs[baseFaceID] = faceI;
                        }
                        else
                        {
                            // Store the shadow values in the local lists
                            shadSubMeshIDs[faceI] =
                                baseMeshShadSubMeshIDs[baseFaceID];
                            shadPatchIDs[faceI] =
                                baseMeshShadPatchIDs[baseFaceID];
                            shadFaceIDs[faceI] =
                                baseMeshShadFaceIDs[baseFaceID];

                            // Update the shadow lists with the local values
                            const label shadSubMeshID =
                                baseMeshShadSubMeshIDs[baseFaceID];
                            //const label shadPatchID =
                            //    baseMeshShadPatchID[baseFaceID];
                            const label shadFaceID =
                                baseMeshShadFaceIDs[baseFaceID];

                            interfaceShadowSubMeshID_
                            [
                                shadSubMeshID
                            ][shadFaceID] = matI;
                            interfaceShadowPatchID_[shadSubMeshID][shadFaceID] =
                                patchI;
                            interfaceShadowFaceID_[shadSubMeshID][shadFaceID] =
                                faceI;
                        }
                    }
                    else
                    {
                        // Shadow IDs not set faces that are on a processor
                        // patch in the base mesh
                    }
                }
            }
        }
    }
}


const Foam::PtrList<Foam::labelList>&
Foam::solidSubMeshes::interfaceShadowSubMeshID() const
{
    if (interfaceShadowSubMeshID_.empty())
    {
        calcInterfaceShadowIDs();
    }

    return interfaceShadowSubMeshID_;
}


const Foam::PtrList<Foam::labelList>&
Foam::solidSubMeshes::interfaceShadowPatchID() const
{
    if (interfaceShadowPatchID_.empty())
    {
        calcInterfaceShadowIDs();
    }

    return interfaceShadowPatchID_;
}


const Foam::PtrList<Foam::labelList>&
Foam::solidSubMeshes::interfaceShadowFaceID() const
{
    if (interfaceShadowFaceID_.empty())
    {
        calcInterfaceShadowIDs();
    }

    return interfaceShadowFaceID_;
}


void Foam::solidSubMeshes::makeInterfaceBaseFaces() const
{
    if (interfaceBaseFacesPtr_)
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    if (cellZoneNames_.size() > 1)
    {
        // Create material index field from cellZones

        volScalarField materials
        (
            IOobject
            (
                "materials",
                baseMesh().time().timeName(),
                baseMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            baseMesh(),
            dimensionedScalar("0", dimless, 0)
        );

        scalarField& materialsI = materials.primitiveFieldRef();

        const meshCellZones& cellZones = baseMesh().cellZones();
        forAll(cellZoneNames_, matI)
        {
            const label cellZoneID =
                cellZones.findZoneID(cellZoneNames_[matI]);

            if (cellZoneID < 0)
            {
                FatalErrorInFunction
                    << "cellZone " << cellZoneNames_[matI]
                    << " not found for material " << cellZoneNames_[matI]
                    << abort(FatalError);
            }

            const labelList& curCellZone =
                baseMesh().cellZones()[cellZoneID];

            forAll(curCellZone, cI)
            {
                const label cellID = curCellZone[cI];
                materialsI[cellID] = matI;
            }
        }

        // Sync coupled boundaries
        materials.correctBoundaryConditions();

        const unallocLabelList& owner = baseMesh().owner();
        const unallocLabelList& neighbour = baseMesh().neighbour();

        labelHashSet interFacesSet;

        forAll(neighbour, faceI)
        {
            if
            (
                mag(materialsI[neighbour[faceI]] - materialsI[owner[faceI]])
              > SMALL
            )
            {
                interFacesSet.insert(faceI);
            }
        }

        forAll(materials.boundaryField(), patchI)
        {
            if
            (
                baseMesh().boundary()[patchI].type()
             == processorFvPatch::typeName
            )
            {
                const scalarField ownMat
                (
                    materials.boundaryField()[patchI].patchInternalField()
                );

                const scalarField ngbMat
                (
                    materials.boundaryField()[patchI].patchNeighbourField()
                );

                forAll(ownMat, faceI)
                {
                    if (mag(ownMat[faceI] - ngbMat[faceI]) > SMALL)
                    {
                        const label globalFaceID =
                            baseMesh().boundaryMesh()[patchI].start() + faceI;

                        interFacesSet.insert(globalFaceID);
                    }
                }
            }
        }

        interfaceBaseFacesPtr_ = new labelList(interFacesSet.toc());
    }
    else
    {
        interfaceBaseFacesPtr_ = new labelList(0);
    }
}


const Foam::labelList& Foam::solidSubMeshes::interfaceBaseFaces() const
{
    if (!interfaceBaseFacesPtr_)
    {
        makeInterfaceBaseFaces();
    }

    return *interfaceBaseFacesPtr_;
}


void Foam::solidSubMeshes::makePointNumOfMaterials() const
{
    if (pointNumOfMaterialsPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    pointNumOfMaterialsPtr_ = new labelList(baseMesh().nPoints(), 0);
    labelList& pointNumOfMaterials = *pointNumOfMaterialsPtr_;

    // Create material index field from cellZones

    scalarField materialsI(baseMesh().nCells(), 0);

    const meshCellZones& cellZones = baseMesh().cellZones();
    forAll(cellZoneNames_, matI)
    {
        const label cellZoneID =
            cellZones.findIndex(cellZoneNames_[matI]);

        if (cellZoneID < 0)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneNames_[matI]
                << " not found for material " << cellZoneNames_[matI]
                << abort(FatalError);
        }

        const labelList& curCellZone =
            baseMesh().cellZones()[cellZoneID];

        forAll(curCellZone, cI)
        {
            const label cellID = curCellZone[cI];
            materialsI[cellID] = matI;
        }
    }

    const labelListList& pointCells = baseMesh().pointCells();

    forAll(pointNumOfMaterials, pointI)
    {
        // Count the number of unique materials in adjacent cells
        const labelList& curCells = pointCells[pointI];

        labelHashSet matSet;

        forAll(curCells, cellI)
        {
            if (!matSet.found(materialsI[curCells[cellI]]))
            {
                matSet.insert(materialsI[curCells[cellI]]);
            }
        }

        pointNumOfMaterials[pointI] = matSet.toc().size();
    }
}


const Foam::labelList& Foam::solidSubMeshes::pointNumOfMaterials() const
{
    if (!pointNumOfMaterialsPtr_)
    {
        makePointNumOfMaterials();
    }

    return *pointNumOfMaterialsPtr_;
}


void Foam::solidSubMeshes::makeIsolatedInterfacePoints() const
{
    if (isolatedInterfacePointsPtr_)
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    // Create material index field from cellZones

    scalarField materialsI(baseMesh().nCells(), 0);

    const meshCellZones& cellZones = baseMesh().cellZones();
    forAll(cellZoneNames_, matI)
    {
        const label cellZoneID =
            cellZones.findIndex(cellZoneNames_[matI]);

        if (cellZoneID < 0)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneNames_[matI]
                << " not found for material " << cellZoneNames_[matI]
                << abort(FatalError);
        }

        const labelList& curCellZone =
            baseMesh().cellZones()[cellZoneID];

        forAll(curCellZone, cI)
        {
            const label cellID = curCellZone[cI];
            materialsI[cellID] = matI;
        }
    }

    pointMesh pMesh(baseMesh());

    pointScalarField pointMaterials
    (
        IOobject
        (
            "pointMaterials",
            baseMesh().time().timeName(),
            baseMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh,
        dimensionedScalar("0", dimless, 0)
    );
    scalarField& pointMaterialsI = pointMaterials.primitiveFieldRef();

    const labelListList& pointCells = baseMesh().pointCells();

    forAll(pointMaterialsI, pointI)
    {
        const labelList& curPointCells = pointCells[pointI];

        forAll(curPointCells, cellI)
        {
            label curCell = curPointCells[cellI];
            pointMaterialsI[pointI] += materialsI[curCell];
        }

        pointMaterialsI[pointI] /= curPointCells.size() + SMALL;
    }

    pointMaterials.correctBoundaryConditions();

    scalarField matInter(pointMaterialsI.size(), 0);

    forAll(pointMaterialsI, pointI)
    {
        const labelList& curPointCells = pointCells[pointI];

        if
        (
            mag(pointMaterialsI[pointI] - materialsI[curPointCells[0]])
          > SMALL
        )
        {
            matInter[pointI] = 1;
        }
    }

    labelHashSet isolatedPointsSet;

    const labelList& noMat = pointNumOfMaterials();

    const labelList& spLabels =
        baseMesh().globalData().sharedPointLabels();

    const labelListList& pointFaces = baseMesh().pointFaces();
    forAll(matInter, pointI)
    {
        if (matInter[pointI] && (noMat[pointI] == 1))
        {
            const bool sharedPoint(findIndex(spLabels, pointI) != -1);

            if (!sharedPoint)
            {
                bool hasProcessorFace = false;

                const labelList& curPointFaces = pointFaces[pointI];
                forAll(curPointFaces, faceI)
                {
                    label faceID = curPointFaces[faceI];
                    label patchID =
                        baseMesh().boundaryMesh().whichPatch(faceID);

                    if (patchID != -1)
                    {
                        if
                        (
                            isA<processorPolyPatch>
                            (
                                baseMesh().boundaryMesh()[patchID]
                            )
                        )
                        {
                            if (findIndex(interfaceBaseFaces(), faceID) == -1)
                            {
                                hasProcessorFace = true;
                                break;
                            }
                        }
                    }
                }

                if (hasProcessorFace)
                {
                    isolatedPointsSet.insert(pointI);
                }
            }
        }
    }

    isolatedInterfacePointsPtr_ = new labelList(isolatedPointsSet.toc());
}


const Foam::labelList& Foam::solidSubMeshes::isolatedInterfacePoints() const
{
    if (!isolatedInterfacePointsPtr_)
    {
        makeIsolatedInterfacePoints();
    }

    return *isolatedInterfacePointsPtr_;
}


const Foam::PtrList<Foam::symmTensorField>&
Foam::solidSubMeshes::interfaceShadowSigma() const
{
    if (interfaceShadowSigma_.empty())
    {
        makeInterfaceShadowSigma();
    }

    return interfaceShadowSigma_;
}


void Foam::solidSubMeshes::makeInterfaceShadowSigma() const
{
    if (!interfaceShadowSigma_.empty())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    interfaceShadowSigma_.setSize(subMeshes().size());

    // Set values for each subMesh
    forAll(subMeshes(), subMeshI)
    {
        const meshSubset& subsetMesh = subMeshes()[subMeshI];
        const fvMesh& subMesh = subsetMesh.subMesh();
        const labelList& patchMap = subsetMesh.patchMap();

        // Find the interface patch for the current subMesh
        // we should store this!

        label patchID = -1;

        forAll(subMesh.boundaryMesh(), pI)
        {
            if (patchMap[pI] == -1)
            {
                patchID = pI;
                break;
            }
        }

        if (patchID == -1)
        {
            // Note: if patchID is still -1, it means that this subMesh does not
            // have any faces on an bi-material interface
            interfaceShadowSigma_.set(subMeshI, new symmTensorField(0));
        }
        else
        {
            interfaceShadowSigma_.set
            (
                subMeshI,
                new symmTensorField
                (
                    subMesh.boundaryMesh()[patchID].size(),
                    symmTensor::zero
                )
            );
        }
    }
}


void Foam::solidSubMeshes::updateInterfaceShadowSigma
(
    const bool useVolFieldSigma
)
{
    if (interfaceShadowSigma_.empty())
    {
        makeInterfaceShadowSigma();
    }

    // Field used for syncing the processor patch values
    volSymmTensorField baseSigmaForSyncing
    (
        IOobject
        (
            "baseSigmaForSyncing",
            baseMesh().time().timeName(),
            baseMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        baseMesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    // Set values for each subMesh
    forAll(subMeshes(), subMeshI)
    {
        const meshSubset& subsetMesh = subMeshes()[subMeshI];
        const fvMesh& subMesh = subsetMesh.subMesh();
        const labelList& patchMap = subsetMesh.patchMap();
        const labelList& faceMap = subsetMesh.faceMap();

        // Find the interface patch for the current subMesh
        // we should store this!

        label patchID = -1;

        forAll(subMesh.boundaryMesh(), pI)
        {
            if (patchMap[pI] == -1)
            {
                patchID = pI;
                break;
            }
        }

        if (patchID == -1)
        {
            // This sub mesh has not faces on a bi-material interface
            continue;
        }

        symmTensorField& resultSigma = interfaceShadowSigma_[subMeshI];

        const labelList& interfaceShadowSubMeshID =
            this->interfaceShadowSubMeshID()[subMeshI];
        const labelList& interfaceShadowPatchID =
            this->interfaceShadowPatchID()[subMeshI];
        const labelList& interfaceShadowFaceID =
            this->interfaceShadowFaceID()[subMeshI];

        // Assemble the shadow stress for each face
        forAll(resultSigma, faceI)
        {
            // Check if the face is not on a processor
            if (interfaceShadowSubMeshID[faceI] != -1)
            {
                // ID of the subMesh on the other side of interface
                const label shadowSubMeshID = interfaceShadowSubMeshID[faceI];

                // ID of the subMesh patch on the other side of interface
                const label shadowPatchID = interfaceShadowPatchID[faceI];

                // ID of the subMesh patch face on the other side of interface
                const label shadowFaceID = interfaceShadowFaceID[faceI];

                // Stress calculated at the other side of the interface
                if (useVolFieldSigma)
                {
                    resultSigma[faceI] =
                        subMeshSigma()
                        [
                            shadowSubMeshID
                        ].boundaryField()[shadowPatchID][shadowFaceID];
                }
                else
                {
                    resultSigma[faceI] =
                        subMeshSigmaf()
                        [
                            shadowSubMeshID
                        ].boundaryField()[shadowPatchID][shadowFaceID];
                }
            }
            else
            {
                // Base face is on a processor boundary

                // Local patch start
                const label start = subMesh.boundaryMesh()[patchID].start();

                // Base mesh face ID
                const label baseFaceID = faceMap[start + faceI];

                // Base mesh patch ID
                const label basePatchID =
                    baseMesh().boundaryMesh().whichPatch(baseFaceID);

                // Base patch start
                const label basePatchStart =
                    baseMesh().boundaryMesh()[basePatchID].start();

                // Base mesh patch local face ID
                const label baseLocalFaceID = baseFaceID - basePatchStart;

                // Base mesh patch faceCells
                const unallocLabelList& faceCells =
                    baseMesh().boundaryMesh()[basePatchID].faceCells();

                // Store local stress on the baseMesh proc patch in the patch
                // internal field
                if (useVolFieldSigma)
                {
                    baseSigmaForSyncing[faceCells[baseLocalFaceID]]
                      = subMeshSigma()
                        [
                            subMeshI
                        ].boundaryField()[patchID][faceI];
                }
                else
                {
                    baseSigmaForSyncing[faceCells[baseLocalFaceID]]
                      = subMeshSigmaf()
                        [
                            subMeshI
                        ].boundaryField()[patchID][faceI];
                }
            }
        }
    }

    // Sync base mesh processor patches
    // This will pass the patch internal field and store it on the neighbour
    // patch
    baseSigmaForSyncing.correctBoundaryConditions();

    // Assemble processor values that have been synced
    forAll(subMeshes(), subMeshI)
    {
        const meshSubset& subsetMesh = subMeshes()[subMeshI];
        const fvMesh& subMesh = subsetMesh.subMesh();
        const labelList& patchMap = subsetMesh.patchMap();
        const labelList& faceMap = subsetMesh.faceMap();

        // Find the interface patch for the current subMesh
        // we should store this!

        label patchID = -1;

        forAll(subMesh.boundaryMesh(), pI)
        {
            if (patchMap[pI] == -1)
            {
                patchID = pI;
                break;
            }
        }

        if (patchID == -1)
        {
            // This sub mesh has not faces on a bi-material interface
            continue;
        }

        symmTensorField& resultSigma = interfaceShadowSigma_[subMeshI];

        const labelList& interfaceShadowSubMeshID =
            this->interfaceShadowSubMeshID()[subMeshI];

        forAll(resultSigma, faceI)
        {
            // Check if the face is on a processor
            if (interfaceShadowSubMeshID[faceI] == -1)
            {
                // The base mesh field will now have the patchNeighbourField
                // values stored on the patch

                // Local patch start
                const label start = subMesh.boundaryMesh()[patchID].start();

                // Base mesh face ID
                const label baseFaceID = faceMap[start + faceI];

                // Base mesh patch ID
                const label basePatchID =
                    baseMesh().boundaryMesh().whichPatch(baseFaceID);

                // Base patch start
                const label basePatchStart =
                    baseMesh().boundaryMesh()[basePatchID].start();

                // Base mesh patch local face ID
                const label baseLocalFaceID = baseFaceID - basePatchStart;

                // Copy patch neighbour field values into the result field
                resultSigma[faceI] =
                    baseSigmaForSyncing.boundaryField()
                    [
                        basePatchID
                    ][baseLocalFaceID];
            }
        }
    }
}


bool Foam::solidSubMeshes::biMaterialInterfaceActive() const
{
    if (!biMaterialInterfaceActivePtr_)
    {
        calcBiMaterialInterfaceActive();
    }

    return *biMaterialInterfaceActivePtr_;
}


void Foam::solidSubMeshes::calcBiMaterialInterfaceActive() const
{
    if (biMaterialInterfaceActivePtr_)
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    biMaterialInterfaceActivePtr_ =
        new bool(returnReduce(interfaceBaseFaces().size(), maxOp<int>()));
}


void Foam::solidSubMeshes::clearOut()
{
    subMeshVolToPoint_.clear();
    subMeshSigma_.clear();
    subMeshSigmaf_.clear();
    subMeshD_.clear();
    subMeshGradD_.clear();
    subMeshGradDf_.clear();
    subMeshPointD_.clear();
    deleteDemandDrivenData(biMaterialInterfaceActivePtr_);
    deleteDemandDrivenData(interfaceBaseFacesPtr_);
    interfaceShadowSubMeshID_.clear();
    interfaceShadowPatchID_.clear();
    interfaceShadowFaceID_.clear();
    interfaceShadowSigma_.clear();
    deleteDemandDrivenData(pointNumOfMaterialsPtr_);
    deleteDemandDrivenData(isolatedInterfacePointsPtr_);

    // Make sure to clear the subMeshes after (not before) clearing the subMesh
    // fields
    subsetMeshes_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidSubMeshes::solidSubMeshes
(
    const fvMesh& baseMesh,
    const wordList& cellZoneNames,
    const bool incremental,
    const bool writeSubMeshes
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            baseMesh.time().timeName(),
            baseMesh,
            IOobject::NO_READ,
            writeSubMeshes ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        )
    ),
    baseMesh_(baseMesh),
    cellZoneNames_(cellZoneNames),
    incremental_(incremental),
    subsetMeshes_(),
    subMeshVolToPoint_(),
    subMeshSigma_(),
    subMeshSigmaf_(),
    subMeshD_(),
    subMeshGradD_(),
    subMeshGradDf_(),
    subMeshPointD_(),
    biMaterialInterfaceActivePtr_(NULL),
    interfaceBaseFacesPtr_(NULL),
    interfaceShadowSubMeshID_(),
    interfaceShadowPatchID_(),
    interfaceShadowFaceID_(),
    interfaceShadowSigma_(),
    pointNumOfMaterialsPtr_(NULL),
    isolatedInterfacePointsPtr_(NULL)
{
    // Construct the sub-meshes
    this->subMeshes();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidSubMeshes::~solidSubMeshes()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::solidSubMeshes::baseMesh() const
{
    return baseMesh_;
}


const Foam::solidSubMeshes::meshSubsetList&
Foam::solidSubMeshes::subMeshes() const
{
    if (subsetMeshes_.empty())
    {
        makeSubMeshes();
    }

    return subsetMeshes_;
}


Foam::solidSubMeshes::meshSubsetList& Foam::solidSubMeshes::subMeshes()
{
    if (subsetMeshes_.empty())
    {
        makeSubMeshes();
    }

    return subsetMeshes_;
}


const Foam::solidSubMeshes::volPointInterpList&
Foam::solidSubMeshes::subMeshVolToPoint() const
{
    if (subMeshVolToPoint_.empty())
    {
        makeSubMeshVolToPoint();
    }

    return subMeshVolToPoint_;
}


Foam::PtrList<Foam::volSymmTensorField>& Foam::solidSubMeshes::subMeshSigma()
{
    if (subMeshSigma_.empty())
    {
        calcSubMeshSigma();
    }

    return subMeshSigma_;
}


const Foam::PtrList<Foam::volSymmTensorField>&
Foam::solidSubMeshes::subMeshSigma() const
{
    if (subMeshSigma_.empty())
    {
        calcSubMeshSigma();
    }

    return subMeshSigma_;
}


Foam::PtrList<Foam::surfaceSymmTensorField>&
Foam::solidSubMeshes::subMeshSigmaf()
{
    if (subMeshSigmaf_.empty())
    {
        calcSubMeshSigmaf();
    }

    return subMeshSigmaf_;
}


const Foam::PtrList<Foam::surfaceSymmTensorField>&
Foam::solidSubMeshes::subMeshSigmaf() const
{
    if (subMeshSigmaf_.empty())
    {
        calcSubMeshSigmaf();
    }

    return subMeshSigmaf_;
}


Foam::PtrList<Foam::volVectorField>& Foam::solidSubMeshes::subMeshD()
{
    if (subMeshD_.empty())
    {
        calcSubMeshD();
    }

    return subMeshD_;
}


const Foam::PtrList<Foam::volVectorField>&
Foam::solidSubMeshes::subMeshD() const
{
    if (subMeshD_.empty())
    {
        calcSubMeshD();
    }

    return subMeshD_;
}


Foam::PtrList<Foam::volTensorField>& Foam::solidSubMeshes::subMeshGradD()
{
    if (subMeshGradD_.empty())
    {
        calcSubMeshGradD();
    }

    return subMeshGradD_;
}


const Foam::PtrList<Foam::volTensorField>&
Foam::solidSubMeshes::subMeshGradD() const
{
    if (subMeshGradD_.empty())
    {
        calcSubMeshGradD();
    }

    return subMeshGradD_;
}


Foam::PtrList<Foam::surfaceTensorField>& Foam::solidSubMeshes::subMeshGradDf()
{
    if (subMeshGradDf_.empty())
    {
        calcSubMeshGradDf();
    }

    return subMeshGradDf_;
}


const Foam::PtrList<Foam::surfaceTensorField>&
Foam::solidSubMeshes::subMeshGradDf() const
{
    if (subMeshGradDf_.empty())
    {
        calcSubMeshGradDf();
    }

    return subMeshGradDf_;
}


Foam::PtrList<Foam::pointVectorField>& Foam::solidSubMeshes::subMeshPointD()
{
    if (subMeshPointD_.empty())
    {
        calcSubMeshPointD();
    }

    return subMeshPointD_;
}


const Foam::PtrList<Foam::pointVectorField>&
Foam::solidSubMeshes::subMeshPointD() const
{
    if (subMeshPointD_.empty())
    {
        calcSubMeshPointD();
    }

    return subMeshPointD_;
}


void Foam::solidSubMeshes::interpolateDtoSubMeshD
(
    const volVectorField& D,
    const bool useVolFieldSigma
)
{
    const meshSubsetList& subMeshes = this->subMeshes();
    if (!biMaterialInterfaceActive())
    {
        // No need for any corrections if there are no bi-material interfaces
        forAll(subMeshes, matI)
        {
            subMeshD()[matI] = subMeshes[matI].interpolate(D);
        }

        return;
    }

    const fvMesh& mesh = this->baseMesh();

    // Update the interface shadow sigma fields
    // This can require parallel communication
    updateInterfaceShadowSigma(useVolFieldSigma);

    forAll(subMeshes, matI)
    {
        volVectorField& subMeshD = this->subMeshD()[matI];

        if (!mesh.template foundObject<volScalarField>("impK"))
        {
            subMeshD = subMeshes[matI].interpolate(D);
            continue;
        }

        const fvMesh& subMesh = subMeshes[matI].subMesh();

        const labelList& faceMap = subMeshes[matI].faceMap();
        const labelList& patchMap = subMeshes[matI].patchMap();
        const labelList& cellMap = subMeshes[matI].cellMap();

        // Store interface field as it is overwritten with the interpolated
        // value by the interpolate function
        vectorField DinterfacePrev(0);
        forAll(subMeshD.boundaryField(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                DinterfacePrev = subMeshD.boundaryField()[patchI];
            }
        }

        // Map the base displacement field to the subMesh; this overwrites
        // the interface with the interpolated values
        subMeshD = subMeshes[matI].interpolate(D);

        // Check if a large strain procedure is being used, if so we must
        // calculate the deformed normals
        // If the deformation gradient field 'F' is found, we will assume it is
        // a large/finite strain procedure
        // For finite strain procedures, we will look up the deformation
        // gradient: relative deformation gradient for updated Lagrangian
        // approaches and the total deformation gradient for total
        // approaches
        // What about uns approaches? I may need to include surfaceField options
        // here
        const bool useDeformedNormals = mesh.foundObject<volTensorField>("F");
        tmp<volTensorField> tFinv;
        tmp<volScalarField> tJ;
        if (useDeformedNormals)
        {
            if (mesh.foundObject<volTensorField>("relF"))
            {
                // Updated Lagrangian approach: use the inverse of the relative
                // deformation gradient
                tFinv = inv(mesh.lookupObject<volTensorField>("relF"));
                tJ = det(tFinv());
            }
            else
            {
                // Total Lagrangian approach: use the inverse of the total
                // deformation gradient
                tFinv = inv(mesh.lookupObject<volTensorField>("F"));
                tJ = det(tFinv);
            }
        }

        // Calculate the new correction to the interface valaues
        forAll(subMeshD.boundaryField(), patchI)
        {
            if (patchMap[patchI] == -1)
            {
                // Interface displacement
                vectorField& Dinterface =
                    subMeshD.boundaryFieldRef()[patchI];

                // Patch start face index
                const label start = subMesh.boundaryMesh()[patchI].start();

                // Base mesh owner cells
                const unallocLabelList& baseOwn = mesh.owner();

                // Base mesh neighbour cells
                const unallocLabelList& baseNei = mesh.neighbour();

                // Base mesh face interpolation weights
                const surfaceScalarField& baseWeights = mesh.weights();
                const scalarField& baseWeightsI =
                    baseWeights.internalField();

                // Base mesh cell centres
                const volVectorField& baseC = mesh.C();
                const vectorField& baseCI = baseC.internalField();

                // Base mesh face area vectors
                const surfaceVectorField& baseSf = mesh.Sf();
                const vectorField& baseSfI = baseSf.internalField();

                // Base mesh face area vector magnitudes
                const surfaceScalarField& baseMagSf = mesh.magSf();
                const scalarField& baseMagSfI = baseMagSf.internalField();

                // Implicit stiffness field
                const volScalarField& K =
                    mesh.lookupObject<volScalarField>("impK");
                const scalarField& KI = K.internalField();

                // Interface face centres
                const vectorField& patchCf = subMesh.boundary()[patchI].Cf();

                // Stress in the current subMesh at the interface
                const symmTensorField* sigmaPatchPtr = NULL;
                if (useVolFieldSigma)
                {
                    sigmaPatchPtr =
                        &(subMeshSigma()[matI].boundaryField()[patchI]);
                }
                else
                {
                    sigmaPatchPtr =
                        &(subMeshSigmaf()[matI].boundaryField()[patchI]);
                }
                const symmTensorField& sigmaPatch = *sigmaPatchPtr;

                const labelList& faceCells =
                    subMesh.boundaryMesh()[patchI].faceCells();

                // Assemble the shadow sigma field: this is the stress
                // calculated from the other side of the interface (in the
                // subMesh on the other side)
                const symmTensorField& interfaceShadSigma =
                    interfaceShadowSigma()[matI];

                // Calculate the interface displacements
                forAll(Dinterface, faceI)
                {
                    // Base mesh face index
                    const label baseFaceID = faceMap[start + faceI];

                    if (mesh.isInternalFace(baseFaceID))
                    {
                        // Base mesh owner cell index
                        const label baseOwnID = baseOwn[baseFaceID];

                        // Base mesh neighbour cell index
                        const label baseNeiID = baseNei[baseFaceID];

                        // Base mesh face interpolation weight
                        const scalar baseW = baseWeightsI[baseFaceID];

                        // Interface unit normal (on base mesh); this may be in
                        // the opposite direction to the subMesh normal
                        vector n = baseSfI[baseFaceID]/baseMagSfI[baseFaceID];
                        if (useDeformedNormals)
                        {
                            // Interpolate Finv and J to the face
                            const tensorField& FinvI =
                                tFinv().internalField();
                            const tensor Finv =
                                baseW*FinvI[baseOwn[baseFaceID]]
                              + (1.0 - baseW)*FinvI[baseOwn[baseFaceID]];

                            const scalarField& JI = tJ().internalField();
                            const scalar J =
                                baseW*JI[baseOwn[baseFaceID]]
                              + (1.0 - baseW)*JI[baseOwn[baseFaceID]];

                            // Nanson's formula
                            // Note: for updated Lagrangian approach, F is the
                            // relative deformation gradient, whereas for total
                            // Lagrangian approaches, it is the total
                            // deformation gradient
                            n = J*Finv.T() & n;
                            n /= mag(n);
                        }

                        // In surfaceInterpolation.C the deltaCoeffs are
                        // calculated as:
                        // 1.0/max(unitArea & delta, 0.05*mag(delta));

                        // Normal distance from the interface to the cell-centre
                        // on side-a
                        scalar da =
                            mag(n & (patchCf[faceI] - baseCI[baseOwnID]));

                        // Normal distance from the interface to the cell-centre
                        // on side-b
                        scalar db =
                            mag(n & (baseCI[baseNeiID] - patchCf[faceI]));

                        // Lookup the stiffness either side of the interface
                        scalar Ka = KI[baseOwnID];
                        scalar Kb = KI[baseNeiID];

                        // The base own cell should be side-a
                        if (baseOwnID != cellMap[faceCells[faceI]])
                        {
                            n = -n;
                            Swap(Ka, Kb);
                            Swap(da, db);
                        }

                        // Calculate the traction at side-a and side-b
                        const vector tractiona = n & sigmaPatch[faceI];
                        const vector tractionb = n & interfaceShadSigma[faceI];

                        // Calculate the displacement at the interface
                        Dinterface[faceI] =
                            DinterfacePrev[faceI]
                          + (da*db/(db*Ka + da*Kb))*(tractionb - tractiona);
                    }
                    else
                    {
                        // These are faces that are on a processor patch in the
                        // baseMesh but were placed in the oldInternalFaces
                        // patch in the subMesh: these faces are on a
                        // bi-material interface

                        // nei: sb => stored on the patch
                        // nei: db => also stored on the patch: to check
                        // So if I find that the baseFaceID is on a proc patch,
                        // then I can lookup the patch ID and directly lookup
                        // sb and db: nothing else is required

                        // Base mesh patch ID
                        const label basePatchID =
                            mesh.boundaryMesh().whichPatch(baseFaceID);

                        // Base patch start
                        const label basePatchStart =
                            mesh.boundaryMesh()[basePatchID].start();

                        // Base mesh patch local face ID
                        const label baseLocalFaceID =
                            baseFaceID - basePatchStart;

                        // Base mesh owner cell index
                        //const label baseOwnID = baseOwn[baseFaceID];
                        const label baseOwnID =
                            mesh.boundaryMesh()
                            [
                                basePatchID
                            ].faceCells()[baseLocalFaceID];

                        // Base mesh face interpolation weight
                        //const scalar baseW = baseWeightsI[baseFaceID];
                        const scalar baseW =
                            baseWeights.boundaryField()
                            [
                                basePatchID
                            ][baseLocalFaceID];

                        // Interface unit normal (on base mesh); this may be in
                        // the opposite direction to the subMesh normal
                        //vector n = baseSf[baseFaceID]/baseMagSf[baseFaceID];
                        vector n =
                            baseSf.boundaryField()[basePatchID][baseLocalFaceID]
                           /baseMagSf.boundaryField()
                            [
                                basePatchID
                            ][baseLocalFaceID];

                        if (useDeformedNormals)
                        {
                            // Interpolate Finv and J to the face
                            // to-do
                            const tensorField& FinvI =
                                tFinv().internalField();
                            const tensor Finv =
                                baseW*FinvI[baseOwnID]
                              + (1.0 - baseW)
                               *tFinv().boundaryField()
                                [
                                    basePatchID
                                ][baseLocalFaceID];

                            const scalarField& JI = tJ().internalField();
                            const scalar J =
                                baseW*JI[baseOwnID]
                              + (1.0 - baseW)
                               *tJ().boundaryField()
                                [
                                    basePatchID
                                ][baseLocalFaceID];

                            // Nanson's formula
                            // Note: for updated Lagrangian approach, F is the
                            // relative deformation gradient, whereas for total
                            // Lagrangian approaches, it is the total
                            // deformation gradient
                            n = J*Finv.T() & n;
                            n /= mag(n);
                        }

                        // Normal distance from the interface to the cell-centre
                        // on side-a
                        const scalar da =
                            mag(n & (patchCf[faceI] - baseCI[baseOwnID]));

                        // Normal distance from the interface to the cell-centre
                        // on side-b
                        //const scalar db =
                        //    mag(n & (baseCI[baseNeiID] - Cf[faceI]));
                        // Note: processor patches store the patchNeighbourField
                        // directly on the patch, so the patch value will
                        // correspond to the patchNeighbourField value
                        const scalar db =
                            mag
                            (
                                n
                              & (
                                  baseC.boundaryField()
                                  [
                                      basePatchID
                                  ][baseLocalFaceID]
                                - patchCf[faceI]
                              )
                            );

                        // Lookup the stiffness at either side
                        const scalar Ka = KI[baseOwnID];
                        const scalar Kb =
                            K.boundaryField()[basePatchID][baseLocalFaceID];

                        // Calculate the traction at side-a and side-b
                        const vector tractiona = n & sigmaPatch[faceI];
                        const vector tractionb = n & interfaceShadSigma[faceI];

                        // Calculate the displacement at the interface
                        Dinterface[faceI] =
                            DinterfacePrev[faceI]
                          + (da*db/(db*Ka + da*Kb))*(tractionb - tractiona);
                    }
                }
            }
            else
            {
                // These are other patches including processor patches, but it
                // seems that interface processor patch faces are placed in the
                // oldInternalFaces patch, so we check for them above.
                // For now, no need to do anything here.
            }
        }
    }
}


void Foam::solidSubMeshes::correctBoundarySnGrad
(
    PtrList<volVectorField>& subMeshDList,
    PtrList<volTensorField>& subMeshGradDList
)
{
    const meshSubsetList& subMeshes = this->subMeshes();

    forAll(subMeshes, matI)
    {
        const volVectorField& subMeshD = subMeshDList[matI];
        volTensorField& subMeshGradD = subMeshGradDList[matI];
        const fvMesh& subMesh = subMeshes[matI].subMesh();

        forAll(subMeshGradD.boundaryField(), patchI)
        {
            const polyPatch& ppatch = subMesh.boundaryMesh()[patchI];
            const vectorField n(subMesh.boundary()[patchI].nf());
            tensorField& patchGradD =
                subMeshGradD.boundaryFieldRef()[patchI];
            const tensorField patchGradDif
            (
                subMeshGradD.boundaryField()[patchI].patchInternalField()
            );
            const vectorField& patchD = subMeshD.boundaryField()[patchI];
            const vectorField patchDif
            (
                subMeshD.boundaryField()[patchI].patchInternalField()
            );

            // Calculate the corrected snGrad
            vectorField correctedSnGrad;

            if (ppatch.type() == wedgePolyPatch::typeName)
            {
                // Take a reference to the wedge patch
                const wedgePolyPatch& wedgePatch =
                    refCast<const wedgePolyPatch>(ppatch);

                const vectorField& patchC = ppatch.faceCentres();
                const vector& centreN = wedgePatch.centreNormal();
                const vectorField Cn(subMesh.boundary()[patchI].Cn());
                const scalarField d(((Cn - patchC) & centreN)/(n & centreN));
                const vectorField projC(d*n + patchC);

                // Calculate correction vector which connects actual cell
                // centre to the transformed cell centre
                const vectorField k(projC - Cn);

                const Field<vector> projD(patchDif + (k & patchGradDif));

                // Calculate delta coeffs from proj position on centre plane to
                // transformed projected position
                const scalarField projDeltaCoeff
                (
                    1.0/mag(transform(wedgePatch.cellT(), projC) - projC)
                );

                // Calculate the patch snGrad
                correctedSnGrad =
                    (
                        transform(wedgePatch.cellT(), projD) - projD
                    )*projDeltaCoeff;
            }
            else
            {
                const vectorField delta(subMesh.boundary()[patchI].delta());
                const vectorField k(delta - (sqr(n) & delta));
                const scalarField& deltaCoeffs =
                    subMesh.boundary()[patchI].deltaCoeffs();

                correctedSnGrad =
                    (patchD - (patchDif + (k & patchGradDif)))*deltaCoeffs;
            }

            // Update the patch snGrad
            patchGradD += n*(correctedSnGrad - (n & patchGradD));
        }
    }
}


void Foam::solidSubMeshes::correctBoundarySnGradf
(
    PtrList<volVectorField>& subMeshDList,
    PtrList<surfaceTensorField>& subMeshGradDfList,
    PtrList<volTensorField>& subMeshGradDList
)
{
    const meshSubsetList& subMeshes = this->subMeshes();

    forAll(subMeshes, matI)
    {
        const volVectorField& subMeshD = subMeshDList[matI];
        surfaceTensorField& subMeshGradDf = subMeshGradDfList[matI];
        const volTensorField& subMeshGradD = subMeshGradDList[matI];
        const fvMesh& subMesh = subMeshes[matI].subMesh();

        forAll(subMeshGradDf.boundaryField(), patchI)
        {
            const polyPatch& ppatch = subMesh.boundaryMesh()[patchI];
            const vectorField n(subMesh.boundary()[patchI].nf());
            tensorField& patchGradDf =
                subMeshGradDf.boundaryFieldRef()[patchI];
            const tensorField patchGradDif
            (
                subMeshGradD.boundaryField()[patchI].patchInternalField()
            );
            const vectorField& patchD = subMeshD.boundaryField()[patchI];
            const vectorField patchDif
            (
                subMeshD.boundaryField()[patchI].patchInternalField()
            );

            // Calculate the corrected snGrad
            vectorField correctedSnGrad;

            if (ppatch.type() == wedgePolyPatch::typeName)
            {
                // Take a reference to the wedge patch
                const wedgePolyPatch& wedgePatch =
                    refCast<const wedgePolyPatch>(ppatch);

                const vectorField& patchC = ppatch.faceCentres();
                const vector& centreN(wedgePatch.centreNormal());
                const vectorField Cn(subMesh.boundary()[patchI].Cn());
                const scalarField d(((Cn - patchC) & centreN)/(n & centreN));
                const vectorField projC(d*n + patchC);

                // Calculate correction vector which connects actual cell
                // centre to the transformed cell centre
                const vectorField k(projC - Cn);

                const Field<vector> projD(patchDif + (k & patchGradDif));

                // Calculate delta coeffs from proj position on centre plane to
                // transformed projected position
                const scalarField projDeltaCoeff
                (
                    1.0/mag(transform(wedgePatch.cellT(), projC) - projC)
                );

                // Calculate the patch snGrad
                correctedSnGrad =
                    (
                        transform(wedgePatch.cellT(), projD) - projD
                    )*projDeltaCoeff;
            }
            else
            {
                const vectorField delta(subMesh.boundary()[patchI].delta());
                const vectorField k(delta - (sqr(n) & delta));
                const scalarField& deltaCoeffs =
                    subMesh.boundary()[patchI].deltaCoeffs();

                correctedSnGrad =
                    (patchD - (patchDif + (k & patchGradDif)))*deltaCoeffs;
            }

            // Update the patch snGrad
            patchGradDf += n*(correctedSnGrad - (n & patchGradDf));
        }
    }
}


void Foam::solidSubMeshes::moveSubMeshes()
{
    // Sub-meshes only exist when there is more than one material law
    if (cellZoneNames_.size() > 1)
    {
        forAll(subMeshes(), matI)
        {
            Info<< "    Moving subMesh " << subMeshes()[matI].subMesh().name()
                << endl;

            twoDPointCorrector twoDCorrector(subMeshes()[matI].subMesh());
            pointField newPoints
            (
                subMeshes()[matI].subMesh().points() + subMeshPointD()[matI]
            );
            twoDCorrector.correctPoints(newPoints);

            subMeshes()[matI].subMesh().movePoints(newPoints);
            subMeshes()[matI].subMesh().V00();
            subMeshes()[matI].subMesh().moving(false);
//             subMeshes()[matI].subMesh().changing(false);
            subMeshes()[matI].subMesh().setPhi().writeOpt() =
                IOobject::NO_WRITE;
        }
    }
}

bool Foam::solidSubMeshes::writeData(Ostream& os) const
{
    bool good = true;
    if (cellZoneNames_.size() > 1)
    {
        forAll(subMeshes(), matI)
        {
            good = good && subMeshes()[matI].subMesh().write();
        }
    }
    return good;
}

const Foam::fvMeshSubset& Foam::solidSubMeshes::operator[]
(
    const word& regionName
) const
{
    forAll(subsetMeshes_, regioni)
    {
        if (subsetMeshes_[regioni].subMesh().name() == regionName)
        {
            return subsetMeshes_[regioni];
        }
    }

    wordList regionNames(subsetMeshes_.size());
    forAll(subsetMeshes_, regioni)
    {
        regionNames[regioni] = subsetMeshes_[regioni].subMesh().name();
    }
    FatalErrorInFunction
        << "No region " << regionName << " could be found." << endl
        << "Valid regions are:" << nl
        << regionNames << endl
        << abort(FatalError);
    return subsetMeshes_[0];
}


// ************************************************************************* //
