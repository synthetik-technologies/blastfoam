/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#ifndef immersedBoundaryObjects_H
#define immersedBoundaryObjects_H

#include "immersedBoundaryObject.H"

#include "MovingImmersedBoundaryObject.H"
#include "StationaryImmersedBoundaryObject.H"

#include "TemperatureIndependentImmersedBoundaryObject.H"
#include "IsothermalImmersedBoundaryObject.H"
#include "AnisothermalImmersedBoundaryObject.H"

#include "NondeformableImmersedMesh.H"
#include "DeformableImmersedMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    typedef
        StationaryImmersedBoundaryObject
        <
            TemperatureIndependentImmersedBoundaryObject
            <
                immersedBoundaryObject
            >
        > stationaryImmersedBoundaryObject;

    typedef
        MovingImmersedBoundaryObject
        <
            TemperatureIndependentImmersedBoundaryObject
            <
                immersedBoundaryObject
            >
        > movingImmersedBoundaryObject;

    typedef
        StationaryImmersedBoundaryObject
        <
            IsothermalImmersedBoundaryObject
            <
                immersedBoundaryObject
            >
        > isothermalStationaryImmersedBoundaryObject;

    typedef
        MovingImmersedBoundaryObject
        <
            IsothermalImmersedBoundaryObject
            <
                immersedBoundaryObject
            >
        > isothermalMovingImmersedBoundaryObject;

    typedef
        StationaryImmersedBoundaryObject
        <
            AnisothermalImmersedBoundaryObject
            <
                NondeformableImmersedMesh
                <
                    immersedBoundaryObject
                >
            >
        > anisothermalStationaryImmersedBoundaryObject;

    typedef
        MovingImmersedBoundaryObject
        <
            AnisothermalImmersedBoundaryObject
            <
                NondeformableImmersedMesh
                <
                    immersedBoundaryObject
                >
            >
        > anisothermalMovingImmersedBoundaryObject;

    typedef
        StationaryImmersedBoundaryObject
        <
            AnisothermalImmersedBoundaryObject
            <
                DeformableImmersedMesh
                <
                    immersedBoundaryObject
                >
            >
        > deformableImmersedBoundaryObject;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
