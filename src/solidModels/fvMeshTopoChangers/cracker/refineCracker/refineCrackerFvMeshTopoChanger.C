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

#include "refineCrackerFvMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CrackerMesh>
void Foam::RefineCrackerFvMesh<CrackerMesh>::readDict()
{
    const dictionary refineDict
    (
        this->dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );
    refiner_->readDict(refineDict);
    error_->read(refineDict);
}


template<class CrackerMesh>
void Foam::RefineCrackerFvMesh<CrackerMesh>::distribute
(
    const mapDistributePolyMesh& map
)
{
    refiner_->distribute(map);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CrackerMesh>
Foam::RefineCrackerFvMesh<CrackerMesh>::RefineCrackerFvMesh
(
    const IOobject& io
)
:
    CrackerMesh(io),
    error_(errorEstimator::New(*this, this->dynamicMeshDict())),
    refiner_(fvMeshRefiner::New(*this, this->dynamicMeshDict())),
    curTimeIndex_(-1)
{
    // Add zones and mesh modifiers
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CrackerMesh>
Foam::RefineCrackerFvMesh<CrackerMesh>::~RefineCrackerFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CrackerMesh>
bool Foam::RefineCrackerFvMesh<CrackerMesh>::update()
{
    return refine();
}


template<class CrackerMesh>
bool Foam::RefineCrackerFvMesh<CrackerMesh>::refine()
{
    // Re-read dictionary. Chosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    readDict();

    bool refined = false;
    if (curTimeIndex_ != this->time().timeIndex())
    {
        //- Update error
        error_->update();
        error_->error().correctBoundaryConditions();

        refined =
            refiner_->refine(error_->error(), error_->maxRefinement());
        curTimeIndex_ = this->time().timeIndex();
    }
    bool cracked = CrackerMesh::update();
    return refined || cracked;
}


template<class CrackerMesh>
bool Foam::RefineCrackerFvMesh<CrackerMesh>::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    return
        CrackerMesh::writeObject(fmt, ver, cmp, write)
     && refiner_->write();
}


// ************************************************************************* //
