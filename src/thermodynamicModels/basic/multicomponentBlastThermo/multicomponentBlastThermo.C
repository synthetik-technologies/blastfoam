/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based BasicThermodynamic
                            class
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

#include "multicomponentBlastThermo.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentBlastThermo::multicomponentBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    basicSpecieBlastMixture
    (
        dict,
        dict.lookup<wordList>("species"),
        mesh,
        phaseName
    ),
    mesh_(mesh),
    masterName_(masterName),
    normalise_(masterName != word::null),
    massTransferRates_(this->species_.size()),
    implicitSources_(this->species_.size())
{}


Foam::multicomponentBlastThermo::multicomponentBlastThermo
(
    const speciesTable& species,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    basicSpecieBlastMixture
    (
        dict,
        species,
        mesh,
        phaseName
    ),
    mesh_(mesh),
    masterName_(masterName),
    normalise_(masterName != word::null),
    massTransferRates_(this->species_.size()),
    implicitSources_(this->species_.size())
{}


Foam::multicomponentBlastThermo::integrator::integrator
(
    const fvMesh& mesh,
    PtrList<volScalarField>& Y,
    PtrListDictionary<volScalarField::Internal>& massTransferRates,
    PtrListDictionary<fvScalarMatrix>& implicitSources,
    const List<bool>& active,
    const word& alphaRhoName,
    const word& alphaRhoPhiName,
    const bool normalize
)
:
    timeIntegrationSystem
    (
        IOobject::groupName("multicomponentIntegrator", Y[0].group()),
        mesh
    ),
    mesh_(mesh),
    Y_(Y),
    massTransferRates_(massTransferRates),
    implicitSources_(implicitSources),
    active_(active),
    alphaRho_(mesh_.lookupObject<volScalarField>(alphaRhoName)),
    alphaRhoPhi_(mesh_.lookupObject<surfaceScalarField>(alphaRhoPhiName)),
    normalize_(normalize)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multicomponentBlastThermo::~multicomponentBlastThermo()
{}

Foam::multicomponentBlastThermo::integrator::~integrator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multicomponentBlastThermo::initializeModels()
{
    const fvMesh& mesh = Y_[0].mesh();
    word alphaRhoName(IOobject::groupName("rho", masterName_));
    word alphaRhoPhiName(IOobject::groupName("rhoPhi", masterName_));
    if
    (
        mesh.foundObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", phaseName_)
        )
    )
    {
        alphaRhoPhiName = IOobject::groupName("alphaRhoPhi", phaseName_);
    }
    else if
    (
        mesh.foundObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", masterName_)
        )
    )
    {
        alphaRhoPhiName = IOobject::groupName("alphaRhoPhi", masterName_);
    }

    if
    (
        mesh.foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", phaseName_)
        )
    )
    {
        alphaRhoName = IOobject::groupName("alphaRho", phaseName_);
    }
    else if
    (
        mesh.foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", masterName_)
        )
    )
    {
        alphaRhoName = IOobject::groupName("alphaRho", masterName_);
    }

    integratorPtr_.reset
    (
        new integrator
        (
            mesh,
            Y_,
            massTransferRates_,
            implicitSources_,
            active_,
            alphaRhoName,
            alphaRhoPhiName,
            defaultSpeciei_ < 0
        )
    );
}


void Foam::multicomponentBlastThermo::update()
{
    integratorPtr_->update();
    if (integratorPtr_->normalize())
    {
        defaultSpeciei_ = -1;
    }
    correctMassFractions();
}


void Foam::multicomponentBlastThermo::solve()
{
    integratorPtr_->solve();
    if (integratorPtr_->normalize())
    {
        defaultSpeciei_ = -1;
    }
}


void Foam::multicomponentBlastThermo::postUpdate()
{
    integratorPtr_->postUpdate();
    if (integratorPtr_->normalize())
    {
        defaultSpeciei_ = -1;
    }
    correctMassFractions();
}


void Foam::multicomponentBlastThermo::clearDeltas()
{
    forAll(massTransferRates_, i)
    {
        if (massTransferRates_.PtrList<volScalarField::Internal>::set(i))
        {
            massTransferRates_[i] = Zero;
        }
    }
}


void Foam::multicomponentBlastThermo::clearSources()
{
    forAll(implicitSources_, i)
    {
        if (implicitSources_.PtrList<fvScalarMatrix>::set(i))
        {
            implicitSources_[i] *= 0;
        }
    }
}


void Foam::multicomponentBlastThermo::addDelta
(
    const word& name,
    tmp<volScalarField>&& delta
)
{
    if (this->containsSpecie(name))
    {
        defaultSpeciei_ = -1;
        const label speciei = species_[name];
        if (!massTransferRates_.PtrList<volScalarField::Internal>::set(speciei))
        {
            massTransferRates_.set
            (
                speciei,
                name,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        "massTransfer:" + IOobject::groupName(name, phaseName_),
                        mesh_.time().name(),
                        mesh_
                    ),
                    delta()()
                )
            );
        }
        else
        {
            massTransferRates_[speciei] += delta();
        }
        delta.clear();
    }
}



void Foam::multicomponentBlastThermo::addDelta
(
    const word& name,
    tmp<volScalarField::Internal>&& delta
)
{
    if (this->containsSpecie(name))
    {
        defaultSpeciei_ = -1;
        const label speciei = species_[name];
        if (!massTransferRates_.PtrList<volScalarField::Internal>::set(speciei))
        {
            massTransferRates_.set
            (
                speciei,
                name,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        "massTransfer:" + IOobject::groupName(name, phaseName_),
                        mesh_.time().name(),
                        mesh_
                    ),
                    delta()
                )
            );
        }
        else
        {
            massTransferRates_[speciei] += delta();
        }
        delta.clear();
    }
}


void Foam::multicomponentBlastThermo::addDelta
(
    const word& name,
    const volScalarField::Internal& delta
)
{
    if (this->containsSpecie(name))
    {
        defaultSpeciei_ = -1;
        const label speciei = species_[name];
        if (!massTransferRates_.PtrList<volScalarField::Internal>::set(speciei))
        {
            massTransferRates_.set
            (
                speciei,
                name,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        "massTransfer:" + IOobject::groupName(name, phaseName_),
                        mesh_.time().name(),
                        mesh_
                    ),
                    delta
                )
            );
        }
        else
        {
            massTransferRates_[speciei] += delta;
        }
    }
}


void Foam::multicomponentBlastThermo::addSource
(
    const word& name,
    tmp<fvScalarMatrix>& source
)
{
    if (species_.found(name))
    {
        defaultSpeciei_ = -1;
        if (implicitSources_.found(name))
        {
            implicitSources_[name] += source;
        }
        else
        {
            implicitSources_.set
            (
                species_[name],
                name,
                source
            );
        }
    }
}


void Foam::multicomponentBlastThermo::addSource
(
    const word& name,
    tmp<fvScalarMatrix>&& source
)
{
    if (species_.found(name))
    {
        defaultSpeciei_ = -1;
        if (implicitSources_.found(name))
        {
            implicitSources_[name] += source;
        }
        else
        {
            implicitSources_.set
            (
                species_[name],
                name,
                source
            );
        }
    }
}


void Foam::multicomponentBlastThermo::integrator::update()
{}


void Foam::multicomponentBlastThermo::integrator::solve()
{
    const dimensionedScalar& dT(mesh_.time().deltaT());
    dimensionedScalar residualAlphaRho(dimDensity, 1e-10);

    const volScalarField& alphaRho = alphaRho_;
    tmp<volScalarField> talphaRho0(max(alphaRho_.prevIter(), residualAlphaRho));
    const volScalarField& alphaRho0 = talphaRho0();

    forAll(Y_, i)
    {
        if (active_[i])
        {
            volScalarField& Y = Y_[i];
            volScalarField deltaAlphaRhoY
            (
                fvc::div
                (
                    alphaRhoPhi_,
                    Y,
                    "div(" + alphaRhoPhi_.name() + ",Yi)"
                )
            );
            if (massTransferRates_.PtrList<volScalarField::Internal>::set(i))
            {
                deltaAlphaRhoY.internalFieldRef() -= massTransferRates_[i];
            }

            // External source so not conservative, normalize mass fractions
            if (this->fvTimeInt_->addDeltaSource(Y.name(), deltaAlphaRhoY))
            {
                normalize_ = true;
            }

            // Not conservative, but alphaRho*Yi is
            volScalarField alphaRhoY(alphaRho*Y);
            this->storeAndBlendOld(alphaRhoY);
            this->storeAndBlendDelta(deltaAlphaRhoY);

            Y = (alphaRhoY - dT*deltaAlphaRhoY)/alphaRho0;
            Y.max(0.0);
            Y.correctBoundaryConditions();
        }
    }
}


void Foam::multicomponentBlastThermo::integrator::postUpdate()
{
    dimensionedScalar residualAlphaRho(dimDensity, 1e-10);

    bool isPhase = alphaRho_.group() != word::null;

    UautoPtr<const fluidMulticomponentThermophysicalTransportModel> thermophysicalTransportPtr;
    if
    (
        mesh_.foundObject<fluidMulticomponentThermophysicalTransportModel>
        (
            IOobject::groupName("thermophysicalTransport", alphaRho_.group())
        )
    )
    {
        thermophysicalTransportPtr.set
        (
            &mesh_.lookupObject<fluidMulticomponentThermophysicalTransportModel>
            (
                IOobject::groupName("thermophysicalTransport", alphaRho_.group())
            )
        );
    }

    forAll(Y_, i)
    {
        volScalarField& Yi(Y_[i]);

        bool needUpdate =
            (
                active_[i]
             && (
                    implicitSources_.PtrList<fvScalarMatrix>::set(i)
                 || thermophysicalTransportPtr.valid()
                )
            )
         || this->needSolve(Yi.name());

        if (needUpdate)
        {
            fvScalarMatrix YEqn
            (
                fvm::ddt(alphaRho_, Yi)
              - fvc::ddt(alphaRho_.prevIter(), Yi)
             ==
                models().source(alphaRho_, Yi)
            );
            if (isPhase)
            {
                YEqn +=
                    fvm::ddt(residualAlphaRho, Yi)
                  - fvc::ddt(residualAlphaRho, Yi);
            }

            if (implicitSources_.PtrList<fvScalarMatrix>::set(i))
            {
                YEqn -= implicitSources_[i];
            }

            if (thermophysicalTransportPtr.valid())
            {
                YEqn += thermophysicalTransportPtr->divj(Yi);
            }

            constraints().constrain(YEqn);
            YEqn.solve("Yi");
            Yi.max(0.0);
            constraints().constrain(Yi);

        }
    }
}


// ************************************************************************* //
