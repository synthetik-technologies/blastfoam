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


inline const Foam::word& Foam::immersedMeshMapper::name() const
{
    return immersedObjectPtr_->name();
}


inline const Foam::fvMesh& Foam::immersedMeshMapper::immersedMesh() const
{
    return immersedMesh_;
}


inline const Foam::immersedBoundaryObject&
Foam::immersedMeshMapper::immersedObject() const
{
    return *immersedObjectPtr_;
}


inline Foam::label Foam::immersedMeshMapper::interfaceIndex() const
{
    return interfaceIndex_;
}

// ************************************************************************* //
