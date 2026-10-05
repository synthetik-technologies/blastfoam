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
#include "zeroGradientFvPatchFields.H"

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
    ),
    alphaCutoff_
    (
        dict.lookupOrDefault<scalar>
        (
            "alphaCutoff",
            dict.lookupOrDefault<scalar>("alphaMin", 0.0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::dispersedDragModel::~dispersedDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::dispersedDragModel::Ki() const
{
    tmp<volScalarField> tKi
    (
        new volScalarField
        (
            IOobject
            (
                "Ki",
                pair_.continuous().mesh().time().name(),
                pair_.continuous().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pair_.continuous().mesh(),
            dimensionedScalar(dimDensity/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& ki = tKi.ref();

    forAll(ki, celli)
    {
        ki[celli] = cellKi(celli);
    }
    ki.correctBoundaryConditions();
    return tKi;
}


Foam::scalar Foam::dragModels::dispersedDragModel::cellKi
(
    const label celli
) const
{
    const scalar d = pair_.dispersed().celld(celli);
    const scalar d0 = pair_.dispersed().celld0(celli);
    const scalar dEff = max(d, d0);
    return
        0.75
       *cellCdRe(celli)
       *swarmCorrection_->cellCs(celli)
       *pair_.continuous().rho()[celli]
       *max(pair_.continuous().cellnu(celli), 1e-15)
       *(d*d)
       /pow4(dEff);
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::dispersedDragModel::stabK() const
{
    tmp<volScalarField> tstabK
    (
        new volScalarField
        (
            IOobject
            (
                "stabK",
                pair_.continuous().mesh().time().name(),
                pair_.continuous().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pair_.continuous().mesh(),
            dimensionedScalar(dimDensity/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& stabK = tstabK.ref();

    const volScalarField& alphaD = pair_.dispersed();
    const scalar rAlphaD = pair_.dispersed().residualAlpha().value();
    const volScalarField& alphaC = pair_.continuous();
    const scalar rAlphaC = pair_.continuous().residualAlpha().value();
    const scalar cutoffD = max(rAlphaD, alphaCutoff_);
    const scalar cutoffC = max(rAlphaC, alphaCutoff_);

    forAll(stabK, celli)
    {
        if (alphaD[celli] > cutoffD && alphaC[celli] > cutoffC)
        {
            stabK[celli] = max(alphaD[celli], rAlphaD)*cellKi(celli);
        }
    }
    stabK.correctBoundaryConditions();
    return tstabK;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::dispersedDragModel::stabKf() const
{
    const scalar cutoffD = max(pair_.dispersed().residualAlpha().value(), alphaCutoff_);
    const scalar cutoffC = max(pair_.continuous().residualAlpha().value(), alphaCutoff_);
    if (alphaCutoff_ > 0.0)
    {
        return
            pos0(fvc::interpolate(pair_.dispersed()) - cutoffD)
           *pos0(fvc::interpolate(pair_.continuous()) - cutoffC)
           *max
            (
                fvc::interpolate(pair_.dispersed()),
                pair_.dispersed().residualAlpha()
            )*fvc::interpolate(Ki());
    }
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
    const scalar cutoffD = max(pair_.dispersed().residualAlpha().value(), alphaCutoff_);
    const scalar cutoffC = max(pair_.continuous().residualAlpha().value(), alphaCutoff_);
    if
    (
        pair_.dispersed()[celli] > cutoffD
     && pair_.continuous()[celli] > cutoffC
    )
    {
        return
            max
            (
                pair_.dispersed()[celli],
                pair_.dispersed().residualAlpha().value()
            )*cellKi(celli);
    }
    return 0.0;
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::dispersedDragModel::K() const
{
    if (alphaCutoff_ > 0.0)
    {
        return
            pos0(pair_.dispersed() - alphaCutoff_)
           *pos0(pair_.continuous() - alphaCutoff_)
           *pair_.dispersed()*Ki();
    }
    return pair_.dispersed()*Ki();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::dispersedDragModel::Kf() const
{
    if (alphaCutoff_ > 0.0)
    {
        return
            pos0(fvc::interpolate(pair_.dispersed()) - alphaCutoff_)
           *pos0(fvc::interpolate(pair_.continuous()) - alphaCutoff_)
           *fvc::interpolate(pair_.dispersed())
           *fvc::interpolate(Ki());
    }
    return
        fvc::interpolate(pair_.dispersed())
       *fvc::interpolate(Ki());
}


Foam::scalar Foam::dragModels::dispersedDragModel::cellK
(
    const label celli
) const
{
    const scalar cutoffD = max(pair_.dispersed().residualAlpha().value(), alphaCutoff_);
    const scalar cutoffC = max(pair_.continuous().residualAlpha().value(), alphaCutoff_);
    if (pair_.dispersed()[celli] <= cutoffD || pair_.continuous()[celli] <= cutoffC)
    {
        return 0.0;
    }
    return pair_.dispersed()[celli]*cellKi(celli);
}


// ************************************************************************* //
