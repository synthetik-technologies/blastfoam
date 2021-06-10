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

#include "multicomponentThermoModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo>
Foam::multicomponentThermoModel<BasicThermo>::multicomponentThermoModel
(
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    BasicThermo
    (
        name,
        p,
        rho,
        e,
        T,
        dict,
        master,
        masterName
    ),
    species_(dict.lookup("species")),
    Ys_(species_.size()),
    massTransferRates_(species_.size()),
    inertIndex_
    (
        dict.found("inertSpecie")
      ? species_[dict.lookup<word>("inertSpecie")]
      : -1
    ),
    active_(species_.size(), true)
{
    tmp<volScalarField> tYdefault;
    volScalarField YTot(volScalarField::New("YTot", p.mesh(), 0.0));

    const fvMesh& mesh = rho.mesh();
    forAll(species_, i)
    {
        IOobject header
        (
            IOobject::groupName(species_[i], name),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.typeHeaderOk<volScalarField>(true))
        {
            Ys_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            // Read Ydefault if not already read
            if (!tYdefault.valid())
            {
                word YdefaultName(IOobject::groupName("Ydefault", name));

                IOobject timeIO
                (
                    YdefaultName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                IOobject time0IO
                (
                    YdefaultName,
                    Time::timeName(0),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (timeIO.typeHeaderOk<volScalarField>(true))
                {
                    tYdefault = new volScalarField(timeIO, mesh);
                }
                else
                {
                    tYdefault = new volScalarField(time0IO, mesh);
                }
            }

            Ys_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tYdefault()
                )
            );
        }
        YTot += Ys_[i];

        massTransferRates_.set
        (
            i,
            species_[i],
            new volScalarField
            (
                IOobject
                (
                    "massTransfer:" + IOobject::groupName(species_[i], name),
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", dimDensity/dimTime, 0.0)
            )
        );
    }

    //- Normalize species
    forAll(Ys_, i)
    {
        Ys_[i] /= YTot;
        Ys_[i].correctBoundaryConditions();
    }
}


template<class BasicThermo>
Foam::multicomponentThermoModel<BasicThermo>::multicomponentThermoModel
(
    const speciesTable& species,
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    BasicThermo
    (
        name,
        p,
        rho,
        e,
        T,
        dict,
        master,
        masterName
    ),
    species_(species),
    Ys_(species_.size()),
    massTransferRates_(species_.size()),
    inertIndex_(species_[dict.lookup<word>("inertSpecie")]),
    active_(species_.size(), true)
{
    tmp<volScalarField> tYdefault;
    volScalarField YTot(volScalarField::New("YTot", p.mesh(), 0.0));

    const fvMesh& mesh = rho.mesh();
    forAll(species_, i)
    {
        IOobject header
        (
            IOobject::groupName(species_[i], name),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.typeHeaderOk<volScalarField>(true))
        {
            Ys_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            // Read Ydefault if not already read
            if (!tYdefault.valid())
            {
                word YdefaultName(IOobject::groupName("Ydefault", name));

                IOobject timeIO
                (
                    YdefaultName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                IOobject time0IO
                (
                    YdefaultName,
                    Time::timeName(0),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (timeIO.typeHeaderOk<volScalarField>(true))
                {
                    tYdefault = new volScalarField(timeIO, mesh);
                }
                else
                {
                    tYdefault = new volScalarField(time0IO, mesh);
                }
            }

            Ys_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tYdefault()
                )
            );
        }
        YTot += Ys_[i];
        massTransferRates_.set
        (
            i,
            species_[i],
            new volScalarField
            (
                IOobject
                (
                    "massTransfer:" + IOobject::groupName(species_[i], name),
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", dimDensity/dimTime, 0.0)
            )
        );
    }

    //- Normalize species
    forAll(Ys_, i)
    {
        Ys_[i] /= YTot;
        Ys_[i].correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo>
Foam::multicomponentThermoModel<BasicThermo>::~multicomponentThermoModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo>
void Foam::multicomponentThermoModel<BasicThermo>::initializeModels()
{
    if
    (
        this->rho_.mesh().template foundObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", this->name_)
        )
    )
    {
        alphaRhoPhiName_ = IOobject::groupName("alphaRhoPhi", this->name_);
    }
    else if
    (
        this->rho_.mesh().template foundObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", this->masterName_)
        )
    )
    {
        alphaRhoPhiName_ = IOobject::groupName("alphaRhoPhi", this->masterName_);
    }
    else
    {
        alphaRhoPhiName_ = IOobject::groupName("rhoPhi", this->masterName_);
    }

    if
    (
        this->rho_.mesh().template foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", this->name_)
        )
    )
    {
        alphaRhoName_ = IOobject::groupName("alphaRho", this->name_);
    }
    else if
    (
        this->rho_.mesh().template foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", this->masterName_)
        )
    )
    {
        alphaRhoName_ = IOobject::groupName("alphaRho", this->masterName_);
    }
    else
    {
        alphaRhoName_ = IOobject::groupName("rho", this->masterName_);
    }
}


template<class BasicThermo>
void Foam::multicomponentThermoModel<BasicThermo>::addDelta
(
    const word& name,
    const volScalarField& delta
)
{
    if (massTransferRates_.found(name))
    {
        massTransferRates_[name] += delta;
        return;
    }
    if (this->debug)
    {
        WarningInFunction << name << " was not added to delta." << endl;
    }
}

template<class BasicThermo>
void Foam::multicomponentThermoModel<BasicThermo>::solve()
{
    const dimensionedScalar& dT(this->rho_.mesh().time().deltaT());
    volScalarField YTot(volScalarField::New("YTot", this->rho_.mesh(), 0.0));
    const surfaceScalarField& alphaRhoPhi
    (
        this->rho_.mesh().template lookupObject<surfaceScalarField>
        (
            alphaRhoPhiName_
        )
    );
    const volScalarField& alphaRho
    (
        this->rho_.mesh().template lookupObject<volScalarField>(alphaRhoName_)
    );

    forAll(Ys_, i)
    {
        if (active(i) && i != inertIndex_)
        {
            if (this->step() == 1)
            {
                Ys_[i].storeOldTime();
            }

            volScalarField YOld(Ys_[i]);
            this->storeAndBlendOld(YOld);

            volScalarField deltaAlphaRhoY
            (
                fvc::div(alphaRhoPhi, Ys_[i], "div(Yi)")
              + massTransferRates_[species_[i]]
            );
            this->storeAndBlendDelta(deltaAlphaRhoY);

            Ys_[i] =
                (
                    alphaRho.prevIter()*YOld - dT*deltaAlphaRhoY
                )/max(this->residualRho(), alphaRho);
            Ys_[i].max(0.0);
            Ys_[i].correctBoundaryConditions();

            YTot += Ys_[i];

            // Clear mass transfer after adding
            massTransferRates_[species_[i]] = Zero;
        }
    }

    if (inertIndex_ >= 0)
    {
        Ys_[inertIndex_] = 1.0 - YTot;
        Ys_[inertIndex_].max(0);
    }
    else
    {
        forAll(Ys_, i)
        {
            Ys_[i] /= YTot;
        }
    }
}


// ************************************************************************* //
