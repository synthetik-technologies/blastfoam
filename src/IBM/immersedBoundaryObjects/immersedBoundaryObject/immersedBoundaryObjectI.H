/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline const Foam::word& Foam::immersedBoundaryObject::name() const
{
    return patch_.name();
}


inline const Foam::word& Foam::immersedBoundaryObject::keyword() const
{
    return patch_.name();
}


inline const Foam::polyMesh& Foam::immersedBoundaryObject::pMesh() const
{
    return pMesh_;
}


inline const Foam::polyPatch& Foam::immersedBoundaryObject::patch() const
{
    return patch_;
}


inline Foam::label Foam::immersedBoundaryObject::index() const
{
    return patch_.index();
}


inline const Foam::scalarField&
Foam::immersedBoundaryObject::deltaCoeffs() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return deltaCoeffs_;
}


inline const Foam::labelList&
Foam::immersedBoundaryObject::allInternalCells() const
{
    if (!allInternalCellsPtr_)
    {
        calcMapping();
    }
    return *allInternalCellsPtr_;
}


inline const Foam::labelList&
Foam::immersedBoundaryObject::internalCells() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return *internalCellsPtr_;
}


inline const Foam::labelList&
Foam::immersedBoundaryObject::patchInternalCells() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return patchInternalCells_;
}


inline const Foam::labelList&
Foam::immersedBoundaryObject::patchExternalCells() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return patchExternalCells_;
}


inline const Foam::labelList&
Foam::immersedBoundaryObject::patchCells() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return patchCells_;
}


inline const Foam::labelListList&
Foam::immersedBoundaryObject::patchMap() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return patchMap_;
}


inline const Foam::labelList&
Foam::immersedBoundaryObject::boundaryCells() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return *boundaryCellsPtr_;
}


inline const Foam::labelList&
Foam::immersedBoundaryObject::shellCells() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return *shellCellsPtr_;
}


inline const Foam::labelListList&
Foam::immersedBoundaryObject::interpToCells() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return interpToCells_;
}



inline const Foam::scalarListList&
Foam::immersedBoundaryObject::interpToWeights() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return interpToWeights_;
}


inline Foam::label Foam::immersedBoundaryObject::nFaces() const
{
    if (!shape_.valid())
    {
        initialize();
    }
    return this->size();
}


inline const Foam::vectorField& Foam::immersedBoundaryObject::Sf() const
{
    if (!shape_.valid())
    {
        initialize();
    }
    return this->faceAreas();
}


inline const Foam::scalarField&
Foam::immersedBoundaryObject::magSf() const
{
    if (!magSfPtr_)
    {
        magSfPtr_ = new scalarField(mag(this->Sf()));
    }
    return *magSfPtr_;
}


inline Foam::tmp<Foam::vectorField>
Foam::immersedBoundaryObject::nf() const
{
    return Sf()/magSf();
}


inline const Foam::immersedShape& Foam::immersedBoundaryObject::shape() const
{
    if (!shape_.valid())
    {
        initialize();
    }
    return shape_();
}


inline const Foam::boundBox& Foam::immersedBoundaryObject::bounds() const
{
    if (!shape_.valid())
    {
        initialize();
    }
    return shape_->bounds();
}


inline Foam::tensor
Foam::immersedBoundaryObject::invMomentOfInertia() const
{
    return orientation() & (inv(this->momentOfInertia()) & orientation().T());
}


inline const Foam::vectorField& Foam::immersedBoundaryObject::force() const
{
    return force_;
}


inline Foam::vectorField& Foam::immersedBoundaryObject::force()
{
    return force_;
}


inline const Foam::vector& Foam::immersedBoundaryObject::forceEff() const
{
    return forceEff_;
}


inline Foam::vector& Foam::immersedBoundaryObject::forceEff()
{
    return forceEff_;
}


inline const Foam::vector& Foam::immersedBoundaryObject::momentEff() const
{
    return momentEff_;
}


inline Foam::vector& Foam::immersedBoundaryObject::momentEff()
{
    return momentEff_;
}


inline const Foam::vector& Foam::immersedBoundaryObject::forceExt() const
{
    return forceExt_;
}


inline Foam::vector& Foam::immersedBoundaryObject::forceExt()
{
    return forceExt_;
}


inline const Foam::vector& Foam::immersedBoundaryObject::momentExt() const
{
    return momentExt_;
}


inline Foam::vector& Foam::immersedBoundaryObject::momentExt()
{
    return momentExt_;
}


// ************************************************************************* //
