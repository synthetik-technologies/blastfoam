/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fluxScheme& Foam::compressibleSystem::flux() const
{
    return fluxScheme_();
}


bool Foam::compressibleSystem::inviscid() const
{
    return !turbulence_.valid();
}


const Foam::compressible::momentumTransportModel&
Foam::compressibleSystem::turbulence() const
{
    return turbulence_();
}


Foam::compressible::momentumTransportModel&
Foam::compressibleSystem::turbulence()
{
    return turbulence_();
}


const Foam::fluidThermophysicalTransportModel&
Foam::compressibleSystem::thermophysicalTransport() const
{
    return thermophysicalTransport_();
}


Foam::fluidThermophysicalTransportModel&
Foam::compressibleSystem::thermophysicalTransport()
{
    return thermophysicalTransport_();
}


const Foam::volVectorField& Foam::compressibleSystem::U() const
{
    return U_;
}


const Foam::surfaceScalarField& Foam::compressibleSystem::phi() const
{
    return phi_;
}


const Foam::surfaceScalarField& Foam::compressibleSystem::rhoPhi() const
{
    return rhoPhi_;
}


// ************************************************************************* //
