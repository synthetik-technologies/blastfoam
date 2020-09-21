/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 Alberto Passalacqua
     \\/     M anipulation  |
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

#include "velocityPDFODETransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::velocityPDFODETransportModel::
velocityPDFODETransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const word& support
)
:
    PDFODETransportModel(name, dict, mesh),
    quadrature_(name, mesh, support),
    momentAdvection_
    (
        velocityMomentAdvection::New
        (
            quadrature_.subDict("momentAdvection"),
            quadrature_,
            support
        )
    )
{
    this->initializeODEFields();
    this->lookupAndInitialize();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::velocityPDFODETransportModel::~velocityPDFODETransportModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::velocityPDFODETransportModel::update()
{
    momentAdvection_->update();
}


void Foam::PDFTransportModels::velocityPDFODETransportModel::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    //- Set initial values for stepi
    PtrList<volScalarField> momentsOld(nMoments());

    forAll(momentsOld, mi)
    {
        momentsOld.set
        (
            mi, new volScalarField(quadrature_.moments()[mi])
        );
        this->storeOld(stepi, momentsOld[mi], momentsOld_[mi]);
        this->blendOld(stepi, momentsOld[mi], momentsOld_[mi], ai);
    }

    //- Calculate true deltas for stepi
    PtrList<volScalarField> deltaMoments(nMoments());
    forAll(deltaMoments, mi)
    {
        deltaMoments.set
        (
            mi,
            new volScalarField(momentAdvection_->divMoments()[mi])
        );
        this->storeDelta(stepi, deltaMoments[mi], deltaMoments_[mi]);
        this->blendDelta(stepi, deltaMoments[mi], deltaMoments_[mi], bi);
    }

    dimensionedScalar dT = mesh_.time().deltaT();

    forAll(quadrature_.moments(), mi)
    {
        volScalarField& m = quadrature_.moments()[mi];
        m = momentsOld[mi] - dT*deltaMoments[mi];
        m.correctBoundaryConditions();
    }
    quadrature_.updateQuadrature();
}


void Foam::PDFTransportModels::velocityPDFODETransportModel::postUpdate()
{
    // Solve moment transport equations
    updateImplicitMomentSource();

    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(quadrature_.nMoments());

    // Solve moment transport equations
    forAll(quadrature_.moments(), momenti)
    {
        volVelocityMoment& m = quadrature_.moments()[momenti];

        momentEqns.set
        (
            momenti,
            new fvScalarMatrix
            (
                fvm::ddt(m)
              - fvc::ddt(m)
             ==
                implicitMomentSource(m)
            )
        );
    }

    forAll (momentEqns, mEqni)
    {
        momentEqns[mEqni].relax();
        momentEqns[mEqni].solve();
    }

    quadrature_.updateQuadrature();

    if (solveMomentSources())
    {
        this->explicitMomentSource();
    }
}


// ************************************************************************* //
