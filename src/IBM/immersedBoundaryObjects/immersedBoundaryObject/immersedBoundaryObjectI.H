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

inline const Foam::polyMesh& Foam::immersedBoundaryObject::pMesh() const
{
    return pMesh_;
}


inline bool Foam::immersedBoundaryObject::report() const
{
    return report_;
}


inline const Foam::word& Foam::immersedBoundaryObject::name() const
{
    return shape_->name();
}


inline Foam::scalar Foam::immersedBoundaryObject::mass() const
{
    return mass_;
}


inline Foam::label Foam::immersedBoundaryObject::nFaces() const
{
    return shape_->patch().size();
}


inline const Foam::standAlonePatch&
Foam::immersedBoundaryObject::patch() const
{
    return patch_;
}


inline const Foam::vectorField&
Foam::immersedBoundaryObject::faceCentres() const
{
    return patch_.faceCentres();
}


inline const Foam::pointField& Foam::immersedBoundaryObject::points() const
{
    return patch_.localPoints();
}


inline const Foam::pointField&
Foam::immersedBoundaryObject::faceCentresOld() const
{
    return shape_->faceCentresOld();
}


inline const Foam::vectorField& Foam::immersedBoundaryObject::Sf() const
{
    return patch_.faceAreas();
}


inline Foam::tmp<Foam::scalarField> Foam::immersedBoundaryObject::magSf() const
{
    return mag(Sf());
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


inline const Foam::vectorField&
Foam::immersedBoundaryObject::internalC() const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    return internalC_;
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


inline const Foam::immersedShape& Foam::immersedBoundaryObject::shape() const
{
    return shape_();
}


inline Foam::point Foam::immersedBoundaryObject::initialCentreOfMass() const
{
    return shape_->centreOfMass();
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


inline const Foam::vector& Foam::immersedBoundaryObject::g() const
{
    return g_;
}


inline const Foam::boundBox& Foam::immersedBoundaryObject::bounds() const
{
    return shape_->bounds();
}

// ************************************************************************* //
