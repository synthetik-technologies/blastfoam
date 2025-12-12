/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
2025-06-09 Jeff Heylmun:    Added cell based returns
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

#include "dispersedDragModel.H"
#include "BlendedInterfacialModel.H"
#include "phasePair.H"
#include "noSwarmCorrection.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(dispersedDragModel, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::dragModels::dispersedDragModel::dispersedDragModel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    swarmCorrection_
    (
        dict.isDict(swarmCorrection::typeName)
      ? swarmCorrection::New(dict, pair)
      : autoPtr<swarmCorrection>(new swarmCorrections::noSwarm(dict, pair))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::dispersedDragModel::~dispersedDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::dispersedDragModel::Ki() const
{
    return
        0.75
       *CdRe()
       *swarmCorrection_->Cs()
       *pair_.continuous().rho()
       *pair_.continuous().nu()
       /sqr(pair_.dispersed().d());
}


Foam::scalar Foam::dragModels::dispersedDragModel::cellKi
(
    const label celli
) const
{
    return
        0.75
       *cellCdRe(celli)
       *swarmCorrection_->cellCs(celli)
       *pair_.continuous().rho()[celli]
       *pair_.continuous().cellnu(celli)
       /sqr(pair_.dispersed().celld(celli));
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::dispersedDragModel::stabK() const
{
    return
        max
        (
            pair_.dispersed(),
            pair_.dispersed().residualAlpha()
        )*Ki();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::dispersedDragModel::stabKf() const
{
    return
        max
        (
            fvc::interpolate(pair_.dispersed()),
            pair_.dispersed().residualAlpha()
        )*fvc::interpolate(Ki());
}


Foam::scalar Foam::dragModels::dispersedDragModel::cellStabK
(
    const label celli
) const
{
    return
        max
        (
            pair_.dispersed()[celli],
            pair_.dispersed().residualAlpha().value()
        )*cellKi(celli);
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::dispersedDragModel::K() const
{
    return pair_.dispersed()*Ki();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::dispersedDragModel::Kf() const
{
    return
        fvc::interpolate(pair_.dispersed())
       *fvc::interpolate(Ki());
}


Foam::scalar Foam::dragModels::dispersedDragModel::cellK
(
    const label celli
) const
{
    return pair_.dispersed()[celli]*cellKi(celli);
}


// ************************************************************************* //
