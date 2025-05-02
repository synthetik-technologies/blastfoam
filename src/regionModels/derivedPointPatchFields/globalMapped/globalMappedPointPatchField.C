/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "globalMappedPointPatchField.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "coupledGlobalPolyPatch.H"
#include "vtkWritePolyData.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                p.boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(iF.name()),
    useRefState_(false),
    f0_(),
    refSet_(false)
{
    Field<Type>::operator=(Zero);
}


template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                p.boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(dict.lookup<word>("nbrName")),
    useRefState_(dict.lookupOrDefault<bool>("useRefState", false)),
    f0_(),
    refSet_(false)
{
    if (dict.found("refValue"))
    {
        f0_ = Field<Type>("refValue", dict, this->size());
        refSet_ = true;
    }
}


template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const word& nbrName,
    const bool useRefState
)
:
    fixedValuePointPatchField<Type>(p, iF),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                p.boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(nbrName),
    useRefState_(useRefState),
    f0_(),
    refSet_(false)
{
    Field<Type>::operator=(this->patchInternalField());
}


template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const globalMappedPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                p.boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(ptf.nbrName_),
    useRefState_(ptf.useRefState_),
    f0_(ptf.f0_),
    refSet_(ptf.refSet_)
{
    if (useRefState_ && refSet_)
    {
        f0_ = mapper(ptf.f0_);
    }
}


template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const globalMappedPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                ptf.patch().boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(ptf.nbrName_),
    useRefState_(ptf.useRefState_),
    f0_(ptf.f0_),
    refSet_(ptf.refSet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::globalMappedPointPatchField<Type>::map
(
    const pointPatchField<Type>& ppf,
    const fieldMapper& mapper
)
{
    fixedValuePointPatchField<Type>::map(ppf, mapper);
    if (useRefState_ && refSet_)
    {
        const globalMappedPointPatchField& gmppf =
            dynamicCast<const globalMappedPointPatchField>(ppf);
        mapper(f0_, gmppf.f0_);
    }
}


template<class Type>
void Foam::globalMappedPointPatchField<Type>::reset
(
    const pointPatchField<Type>& ppf
)
{
    fixedValuePointPatchField<Type>::reset(ppf);
    if (useRefState_ && refSet_)
    {
        const globalMappedPointPatchField& gmppf =
            dynamicCast<const globalMappedPointPatchField>(ppf);
        f0_.reset(gmppf.f0_);
    }
}


template<class Type>
void Foam::globalMappedPointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const coupledGlobalPolyPatch& cgpp =
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                this->patch().boundaryMesh().mesh().thisDb()
            )
        )(this->patch());
    const polyMesh& nbrMesh = cgpp.sampleMesh();
    const coupledGlobalPolyPatch& samplePatch = cgpp.samplePatch();
    const label samplePatchi = samplePatch.patch().index();

    const pointPatchField<Type>& pfNbr =
        nbrMesh.lookupObject<GeometricField<Type, pointPatchField, pointMesh>>
        (
            nbrName_
        ).boundaryField()[samplePatchi];
    tmp<Field<Type>> tnbr;
    if (isA<valuePointPatchField<Type>>(pfNbr))
    {
        tnbr = tmp<Field<Type>>
        (
            dynamicCast<const valuePointPatchField<Type>>(pfNbr)
        );
    }
    else
    {
        tnbr = pfNbr.patchInternalField();
    }

//     if (debug > 1 || (debug && this->db().time().outputTime()))
//     {
//         Field<Type> pfGlobal(samplePatch.patchPointToGlobal(nbr));
//         Field<Type> pfInterp
//         (
//             cgpp.patchToPatchInterpolator().transferPoints
//             (
//                 samplePatch.globalPatch(),
//                 pfGlobal
//             )
//         );
//
//         if (Pstream::master())
//         {
//             fileName path
//             (
//                 this->db().time().globalPath()
//                /"VTK"
//                /this->db().time().timeName()
//             );
//             mkDir(path);
//             vtkWritePolyData::write
//             (
//                 path/(this->internalField().name() + "_interpolated.vtk"),
//                 this->internalField().name(),
//                 true,
//                 cgpp.physicalPatch().points(),
//                 labelList(),
//                 edgeList(),
//                 cgpp.physicalPatch(),
//                 this->internalField().name(), true, pfInterp
//             );
//             vtkWritePolyData::write
//             (
//                 path/(nbrName_ + "_actual.vtk"),
//                 nbrName_,
//                 true,
//                 samplePatch.physicalPatch().points(),
//                 labelList(),
//                 edgeList(),
//                 samplePatch.physicalPatch(),
//                 nbrName_,
//                 true,
//                 pfGlobal
//             );
//         }
//     }

    tnbr = samplePatch.pointInterpolate(tnbr);

    if (useRefState_)
    {
        //- Set the reference state if not already set
        if (!refSet_)
        {
            f0_ = tnbr();
            refSet_ = true;
        }

        // Remove reference state
        tnbr.ref() -= f0_;
    }

    Field<Type>::operator=(tnbr);
    fixedValuePointPatchField<Type>::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}


template<class Type>
void Foam::globalMappedPointPatchField<Type>::write(Ostream& os) const
{
    fixedValuePointPatchField<Type>::write(os);
    writeEntry(os, "nbrName", nbrName_);
    writeEntry(os, "useRefState", useRefState_);
    if (useRefState_)
    {
        writeEntry(os, "refValue", f0_);
    }
}

// ************************************************************************* //
