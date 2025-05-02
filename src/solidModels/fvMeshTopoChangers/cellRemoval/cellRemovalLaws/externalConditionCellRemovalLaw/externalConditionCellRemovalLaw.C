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

#include "externalConditionCellRemovalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "globalPolyBoundaryMesh.H"
#include "coupledGlobalPolyPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "removeCells.H"
#include "fvMeshSubset.H"
#include "lookupSolidModel.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(externalConditionCellRemovalLaw, 0);
    addToRunTimeSelectionTable
    (
        cellRemovalLaw, externalConditionCellRemovalLaw, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::externalConditionCellRemovalLaw::externalConditionCellRemovalLaw
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict),
    pCritical_(readScalar(dict.lookup("pCritical"))),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    coupledPatchID_
    (
        mesh.boundaryMesh()[dict.lookup<word>("coupledPatch")].index()
    ),
    patchID_
    (
        mesh.boundaryMesh()[dict.lookup<word>("exposedPatch")].index()
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::externalConditionCellRemovalLaw::~externalConditionCellRemovalLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::externalConditionCellRemovalLaw::cellsToRemove()
{
    const fvPatch& patch = mesh().boundary()[coupledPatchID_];
    const coupledGlobalPolyPatch& cgpp =
        globalPolyBoundaryMesh::New(mesh())(patch.patch());
    const polyMesh& nbrMesh = cgpp.sampleMesh();
    const coupledGlobalPolyPatch& samplePatch = cgpp.samplePatch();
    const label samplePatchi = samplePatch.patch().index();

    const fvPatchField<scalar>& ppNbr =
        nbrMesh.lookupObject<GeometricField<scalar, fvPatchField, volMesh>>
        (
            pName_
        ).boundaryField()[samplePatchi];

    const scalarField pInterp(samplePatch.faceInterpolate(ppNbr));

    const labelList& faceCells = patch.faceCells();
    labelHashSet cellsToRemove;

    forAll(faceCells, fi)
    {
        if (pInterp[fi] > pCritical_)
        {
            cellsToRemove.insert(faceCells[fi]);
        }
    }

    subsetter_.clear();
    meshMap_.clear();
    if (returnReduce(cellsToRemove.size(), sumOp<label>()))
    {
        // Create the subsetter
        subsetter_.set(new fvMeshSubset(mesh()));

        // Subset the mesh
        subsetter_->setLargeCellSubset(cellsToRemove);

        // Rename the mesh
        subsetter_->subMesh().polyMesh::rename(mesh().name() + "_removed");

        // Lookup the solid model
        const solidModel& solid = lookupSolidModel(mesh());

        // If the mesh is not moving, move the subset mesh to the deformed
        // geometry
        if (!solid.movingMesh())
        {
            subsetter_->subMesh().movePoints
            (
                subsetter_->subMesh().points()
                + pointField(solid.pointD(), subsetter_->pointMap())
            );
        }

        #define saveGeoFieldTypes(Type, Patch, Mesh) \
            saveGeoFields<Type, Patch, Mesh>(mesh());

        FOR_ALL_FIELD_TYPES(saveGeoFieldTypes, fvPatchField, volMesh)
        #undef saveGeoFieldTypes
    }

    return cellsToRemove.toc();
}


Foam::label Foam::externalConditionCellRemovalLaw::exposedFacesPatchID()
{
    return patchID_;
}


bool Foam::externalConditionCellRemovalLaw::active() const
{
    return subsetter_.valid();
}


const Foam::meshToMesh& Foam::externalConditionCellRemovalLaw::meshMap() const
{
    if (!meshMap_.valid())
    {
        const fvPatch& patch = mesh().boundary()[coupledPatchID_];
        const coupledGlobalPolyPatch& cgpp =
            globalPolyBoundaryMesh::New(mesh())(patch.patch());
        const polyMesh& nbrMesh = cgpp.sampleMesh();

        //- Create the mesh to mesh mapper using cell volume weighting
        meshMap_.reset
        (
            new meshToMesh
            (
                nbrMesh,
                subsetter_->subMesh(),
                meshToMesh::imCellVolumeWeight,
                false
            )
        );
    }
    return meshMap_();
}


template<>
void Foam::externalConditionCellRemovalLaw::map
(
    const word& fieldName,
    const plusEqOp<scalar>& cop,
    DimensionedField<scalar, volMesh>& field
) const
{
    if (!active())
    {
        return;
    }
    const fvMesh& subMesh = subsetter_->subMesh();
    const fvMesh& mesh = field.mesh();

    const GeometricField<scalar, fvPatchField, volMesh>& rfield =
        subMesh.lookupObject
        <
            GeometricField<scalar, fvPatchField, volMesh>
        >(fieldName);
    const word fieldNameMember(IOobject::member(fieldName));
    const word phaseName(IOobject::group(fieldName));
    if (rfield.dimensions() == dimDensity)
    {
        //- Map mass to new mesh and correct new mass based on cell volume
        const labelListList& addr = meshMap().srcToTgtCellAddr();
        scalar sumV = 0.0;
        forAll(addr, celli)
        {
            if (addr[celli].size())
            {
                sumV += mesh.V()[celli];
            }
        }
        scalar scale(returnReduce(sumV, sumOp<scalar>())/meshMap().V());

        scalarField mass(rfield.primitiveField()*subMesh.V());
        scalarField rhoField(field.size(), Zero);

        meshMap().mapTgtToSrc(mass, cop, rhoField);

        forAll(addr, celli)
        {
            if (addr[celli].size())
            {
                cop(field[celli], rhoField[celli]*scale/mesh.V()[celli]);
            }
        }
    }
    else if
    (
        fieldNameMember == "h" || fieldNameMember == "e"
     || fieldNameMember == "ha" || fieldNameMember == "ea"
    )
    {
        const basicThermo& thermo =
            this->mesh().lookupObject<basicThermo>
            (
                IOobject::groupName(basicThermo::dictName, phaseName)
            );
        const scalarField T(thermo.T(), subsetter_->cellMap());
        const scalarField hs(thermo.hs(T, subsetter_->cellMap()));
        meshMap().mapTgtToSrc(hs, cop, field);
    }
    else
    {
        meshMap().mapTgtToSrc(rfield, cop, field);
    }
}

// ************************************************************************* //
