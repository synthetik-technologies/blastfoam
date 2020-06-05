/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 Alberto Passalacqua
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

#include "univariatePDFODETransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFODETransportModel
::univariatePDFODETransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const surfaceScalarField& phi,
    const word& support
)
:
    PDFODETransportModel(name, dict, mesh),
    quadrature_(name, mesh, support),
    momentAdvection_
    (
        univariateMomentAdvection::New
        (
            quadrature_.subDict("momentAdvection"),
            quadrature_,
            phi,
            support
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFODETransportModel
::~univariatePDFODETransportModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::univariatePDFODETransportModel::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    momentAdvection_().update();

    //- Set initial values for stepi
    PtrList<volScalarField> momentsOld(nMoments());
    forAll(momentsOld, mi)
    {
        momentsOld.set
        (
            mi, new volScalarField(quadrature_.moments()[mi])
        );
    }

    //- Correct of value for changes in cell volume
    if (mesh_.moving() && stepi == 1)
    {
        volScalarField::Internal v0Byv(mesh_.Vsc0()/mesh_.Vsc());
        forAll(momentsOld, mi)
        {
            momentsOld[mi].ref() *= v0Byv;
        }
    }

    // Store old values
    if (oldIs_[stepi - 1] != -1)
    {
        saveFields
        (
            oldIs_[stepi - 1],
            momentsOld,
            momentsOld_
        );
    }

    //- Finish updating initial value for stepi
    forAll(momentsOld, mi)
    {
        momentsOld[mi] *= ai[stepi - 1];
    }
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = oldIs_[i];
        if (fi != -1 && ai[fi] != 0)
        {
            forAll(momentsOld, mi)
            {
                momentsOld[mi] += ai[fi]*momentsOld_[fi][mi];
            }
        }
    }

    //- Store deltas for stepi
    if (deltaIs_[stepi - 1] != -1)
    {
        saveFields
        (
            deltaIs_[stepi - 1],
            momentAdvection_->divMoments(),
            this->deltaMoments_
        );
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
    }
    forAll(deltaMoments, mi)
    {
        deltaMoments[mi] *= bi[stepi - 1];
    }

    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = deltaIs_[i];
        if (fi != -1 && bi[fi] != 0)
        {
            forAll(deltaMoments, mi)
            {
                deltaMoments[mi] += bi[fi]*deltaMoments_[fi][mi];
            }
        }
    }

    dimensionedScalar dT = mesh_.time().deltaT();

    forAll(momentsOld_, mi)
    {
        volScalarField& m = quadrature_.moments()[mi];
        m = momentsOld[mi] - dT*deltaMoments[mi];
        m.correctBoundaryConditions();
    }

    quadrature_.updateQuadrature();

    if (stepi == oldIs_.size())
    {
        // List of moment transport equations
        PtrList<fvScalarMatrix> momentEqns(quadrature_.nMoments());

        // Solve moment transport equations
        forAll(quadrature_.moments(), momenti)
        {
            volScalarMoment& m = quadrature_.moments()[momenti];

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
}


// ************************************************************************* //
