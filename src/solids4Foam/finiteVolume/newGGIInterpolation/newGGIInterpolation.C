/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Description
    Interpolation class dealing with transfer of data between two
    primitivePatches

Author
    Hrvoje Jasak, Wikki Ltd.

Contributor:
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "newGGIInterpolationTemplate.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
label newGGIInterpolation<MasterPatch, SlavePatch>::parMasterStart() const
{
    if (globalData())
    {
        // Integer division intended
        return Foam::min
        (
            masterPatch_.size(),
            Pstream::myProcNo()*(masterPatch_.size()/Pstream::nProcs() + 1)
        );
    }
    else
    {
        // No parallel search: do complete patch
        return 0;
    }
}


template<class MasterPatch, class SlavePatch>
label newGGIInterpolation<MasterPatch, SlavePatch>::parMasterEnd() const
{
    if (globalData())
    {
        // Integer division intended
        return Foam::min
        (
            masterPatch_.size(),
            (Pstream::myProcNo() + 1)*
            (masterPatch_.size()/Pstream::nProcs() + 1)
        );
    }
    else
    {
        // No parallel search: do complete patch
        return masterPatch_.size();
    }
}


template<class MasterPatch, class SlavePatch>
label newGGIInterpolation<MasterPatch, SlavePatch>::parMasterSize() const
{
    return Foam::max
    (
        0,
        this->parMasterEnd() - this->parMasterStart()
    );
}


template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::clearOut()
{
    deleteDemandDrivenData(masterAddrPtr_);
    deleteDemandDrivenData(masterWeightsPtr_);
    deleteDemandDrivenData(slaveAddrPtr_);
    deleteDemandDrivenData(slaveWeightsPtr_);

    deleteDemandDrivenData(uncoveredMasterAddrPtr_);
    deleteDemandDrivenData(uncoveredSlaveAddrPtr_);

    deleteDemandDrivenData(masterPointAddressingPtr_);
    deleteDemandDrivenData(masterPointWeightsPtr_);
    deleteDemandDrivenData(masterPointDistancePtr_);
    deleteDemandDrivenData(masterPointDistanceVectorsPtr_);
    masterEdgeLoopsMap_.clear();

    deleteDemandDrivenData(slavePointAddressingPtr_);
    deleteDemandDrivenData(slavePointWeightsPtr_);
    deleteDemandDrivenData(slavePointDistancePtr_);
    deleteDemandDrivenData(slavePointDistanceVectorsPtr_);
    slaveEdgeLoopsMap_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class MasterPatch, class SlavePatch>
newGGIInterpolation<MasterPatch, SlavePatch>::newGGIInterpolation
(
    const MasterPatch& masterPatch,
    const SlavePatch&  slavePatch,
    const tensorField& forwardT,
    const tensorField& reverseT,
    const vectorField& forwardSep,
    const bool globalData,
    const scalar masterNonOverlapFaceTol,
    const scalar slaveNonOverlapFaceTol,
    const bool rescaleGGIWeightingFactors,
    const quickReject reject,
    const boundBox regionOfInterest
)
:
    masterPatch_(masterPatch),
    slavePatch_(slavePatch),
    forwardT_(forwardT),
    reverseT_(reverseT),
    forwardSep_(forwardSep),
    globalData_(globalData),
    masterNonOverlapFaceTol_(masterNonOverlapFaceTol),
    slaveNonOverlapFaceTol_(slaveNonOverlapFaceTol),
    rescaleGGIWeightingFactors_(rescaleGGIWeightingFactors),
    reject_(reject),
    usePrevCandidateMasterNeighbors_(false),
    prevCandidateMasterNeighbors_(0),
    regionOfInterest_(regionOfInterest),
    masterAddrPtr_(NULL),
    masterWeightsPtr_(NULL),
    masterPointAddressingPtr_(NULL),
    masterPointWeightsPtr_(NULL),
    masterPointDistancePtr_(NULL),
    masterPointDistanceVectorsPtr_(NULL),
    masterEdgeLoopsMap_(0),
    slaveAddrPtr_(NULL),
    slaveWeightsPtr_(NULL),
    slavePointAddressingPtr_(NULL),
    slavePointWeightsPtr_(NULL),
    slavePointDistancePtr_(NULL),
    slavePointDistanceVectorsPtr_(NULL),
    slaveEdgeLoopsMap_(0),
    useNewPointDistanceMethod_(true),
    projectPointsToPatchBoundary_(false),
    checkPointDistanceOrientations_(useNewPointDistanceMethod_ ? false : true),
    uncoveredMasterAddrPtr_(NULL),
    uncoveredSlaveAddrPtr_(NULL)
{
    // Check size of transform.  They should be equal to slave patch size
    // if the transform is not constant
    if (forwardT_.size() > 1 || reverseT_.size() > 1)
    {
        if
        (
            forwardT_.size() != slavePatch_.size()
         || reverseT_.size() != masterPatch_.size()
        )
        {
            FatalErrorIn
            (
                "newGGIInterpolation<MasterPatch, SlavePatch>::"
                "newGGIInterpolation"
            )   << "Incorrectly defined transform: forwardT: "
                << forwardT_.size() << " patch: " << slavePatch_.size()
                << " reverseT: " << reverseT_.size()
                << " patch: " << masterPatch_.size()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
newGGIInterpolation<MasterPatch, SlavePatch>::~newGGIInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const labelListList&
newGGIInterpolation<MasterPatch, SlavePatch>::masterAddr() const
{
    if (!masterAddrPtr_)
    {
        calcAddressing();
    }

    return *masterAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::masterWeights() const
{
    if (!masterWeightsPtr_)
    {
        calcAddressing();
    }

    return *masterWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
newGGIInterpolation<MasterPatch, SlavePatch>::slaveAddr() const
{
    if (!slaveAddrPtr_)
    {
        calcAddressing();
    }

    return *slaveAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::slaveWeights() const
{
    if (!slaveWeightsPtr_)
    {
        calcAddressing();
    }

    return *slaveWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelList&
newGGIInterpolation<MasterPatch, SlavePatch>::uncoveredMasterFaces() const
{
    if (!uncoveredMasterAddrPtr_)
    {
        calcAddressing();
    }

    return *uncoveredMasterAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelList&
newGGIInterpolation<MasterPatch, SlavePatch>::uncoveredSlaveFaces() const
{
    if (!uncoveredSlaveAddrPtr_)
    {
        calcAddressing();
    }

    return *uncoveredSlaveAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
bool newGGIInterpolation<MasterPatch, SlavePatch>::movePoints
(
    const tensorField& forwardT,
    const tensorField& reverseT,
    const vectorField& forwardSep
)
{
    this->forwardT_ = forwardT;
    this->reverseT_ = reverseT;
    this->forwardSep_ = forwardSep;

    if (prevCandidateMasterNeighbors_.size() > 0)
    {
        if (prevCandidateMasterNeighbors_.size() != parMasterSize())
        {
            Info<< "    " << typeName
                << " : clearing prevCandidateMasterNeighbors" << endl;
            clearPrevCandidateMasterNeighbors();
        }
    }

    clearOut();

    return true;
}

template<class MasterPatch, class SlavePatch>
const Foam::List<labelPair>&
newGGIInterpolation<MasterPatch, SlavePatch>::masterPointAddr() const
{
    if (!masterPointAddressingPtr_)
    {
        calcMasterPointAddressing();
    }

    return *masterPointAddressingPtr_;
}

template<class MasterPatch, class SlavePatch>
const Foam::FieldField<Field, scalar>&
newGGIInterpolation<MasterPatch, SlavePatch>::masterPointWeights() const
{
    if (!masterPointWeightsPtr_)
    {
        calcMasterPointWeights();
    }

    return *masterPointWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarField&
newGGIInterpolation<MasterPatch, SlavePatch>
::masterPointDistanceToIntersection() const
{
    if (!masterPointDistancePtr_)
    {
        calcMasterPointAddressing();
    }

    return *masterPointDistancePtr_;
}


template<class MasterPatch, class SlavePatch>
const vectorField&
newGGIInterpolation<MasterPatch, SlavePatch>
::masterPointDistanceVectorsToIntersection() const
{
    if (!masterPointDistanceVectorsPtr_)
    {
        calcMasterPointAddressing();
    }

    return *masterPointDistanceVectorsPtr_;
}


template<class MasterPatch, class SlavePatch>
const Foam::List<labelPair>&
newGGIInterpolation<MasterPatch, SlavePatch>::slavePointAddr() const
{
    if (!slavePointAddressingPtr_)
    {
        calcSlavePointAddressing();
    }

    return *slavePointAddressingPtr_;
}

template<class MasterPatch, class SlavePatch>
const Foam::FieldField<Field, scalar>&
newGGIInterpolation<MasterPatch, SlavePatch>::slavePointWeights() const
{
    if (!slavePointWeightsPtr_)
    {
        calcSlavePointWeights();
    }

    return *slavePointWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarField&
newGGIInterpolation<MasterPatch, SlavePatch>
::slavePointDistanceToIntersection() const
{
    if (!slavePointDistancePtr_)
    {
        calcSlavePointAddressing();
    }

    return *slavePointDistancePtr_;
}


template<class MasterPatch, class SlavePatch>
const vectorField&
newGGIInterpolation<MasterPatch, SlavePatch>
::slavePointDistanceVectorsToIntersection() const
{
    if (!slavePointDistanceVectorsPtr_)
    {
        calcSlavePointAddressing();
    }

    return *slavePointDistanceVectorsPtr_;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> > newGGIInterpolation<MasterPatch, SlavePatch>::
slaveToMasterPointInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != this->slavePatch().nPoints())
    {
        FatalErrorIn
        (
            "newGGIInterpolation::slaveToMasterPointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << this->slavePatch().nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            this->masterPatch().nPoints(),
            pTraits<Type>::zero
        )
    );

    // Escape the interpolation if there are no faces in the target patch
    if (this->masterPatch().nPoints() == 0)
    {
        return tresult;
    }

    Field<Type>& result = tresult();

    const List<typename SlavePatch::FaceType>& slaveFaces =
        this->slavePatch().localFaces();

    const List<labelPair>& addr = masterPointAddr();

    const FieldField<Field, scalar>& weights = masterPointWeights();

    forAll (result, pointI)
    {
        if (addr[pointI].first() > -1)
        {
            const face& hitFace =
                slaveFaces[addr[pointI].first()];

            label pI = addr[pointI].second();

            Type ctrF = average(Field<Type>(pf, hitFace));

            result[pointI] =
                weights[pointI][0]*pf[hitFace[pI]]
              + weights[pointI][1]*pf[hitFace.nextLabel(pI)]
              + weights[pointI][2]*ctrF;
        }
        else
        {
            FatalErrorIn
            (
                "newGGIInterpolation::masterToSlavePointInterpolate"
                "(const Field<Type> pf)"
            )   << "Master point addressing is not correct"
                << abort(FatalError);
        }
    }

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> > newGGIInterpolation<MasterPatch, SlavePatch>::
masterToSlavePointInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != this->masterPatch().nPoints())
    {
        FatalErrorIn
        (
            "newGGIInterpolation::masterToSlavePointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << this->masterPatch().nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            this->slavePatch().nPoints(),
            pTraits<Type>::zero
        )
    );

    // Escape the interpolation if there are no faces in the target patch
    if (this->slavePatch().nPoints() == 0)
    {
        return tresult;
    }

    Field<Type>& result = tresult();

    const List<typename SlavePatch::FaceType>& masterFaces =
        this->masterPatch().localFaces();

    const List<labelPair>& addr = slavePointAddr();

    const FieldField<Field, scalar>& weights = slavePointWeights();

    forAll (result, pointI)
    {
        if (addr[pointI].first() > -1)
        {
            const face& hitFace =
                masterFaces[addr[pointI].first()];

            label pI = addr[pointI].second();

            Type ctrF = average(Field<Type>(pf, hitFace));

            result[pointI] =
                weights[pointI][0]*pf[hitFace[pI]]
              + weights[pointI][1]*pf[hitFace.nextLabel(pI)]
              + weights[pointI][2]*ctrF;
        }
        else
        {
            FatalErrorIn
            (
                "newGGIInterpolation::masterToSlavePointInterpolate"
                "(const Field<Type> pf)"
            )   << "Slave point addressing is not correct"
                << abort(FatalError);
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newGGIInterpolationPolygonIntersection.C"
#   include "newGGIInterpolationQuickRejectTests.C"
#   include "newGGIInterpolationWeights.C"
#   include "newGGIInterpolate.C"

// ************************************************************************* //
