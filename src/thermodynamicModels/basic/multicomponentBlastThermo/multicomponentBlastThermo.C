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
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentBlastThermo::multicomponentBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    integrationSystem
    (
        IOobject::groupName("multicomponentBlastThermo", phaseName),
        mesh
    ),
    mesh_(mesh),
    phaseName_(phaseName),
    masterName_(masterName),
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
    volScalarField YTot(volScalarField::New("YTot", mesh, 0.0));

    forAll(species_, i)
    {
        IOobject header
        (
            IOobject::groupName(species_[i], phaseName),
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
                        IOobject::groupName(species_[i], phaseName),
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
                word YdefaultName(IOobject::groupName("Ydefault", phaseName));

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
                        IOobject::groupName(species_[i], phaseName),
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
                    "massTransfer:" + IOobject::groupName(species_[i], phaseName),
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


Foam::multicomponentBlastThermo::multicomponentBlastThermo
(
    const speciesTable& species,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    integrationSystem
    (
        IOobject::groupName("multicomponentBlastThermo", phaseName),
        mesh
    ),
    mesh_(mesh),
    phaseName_(phaseName),
    masterName_(masterName),
    species_(species),
    Ys_(species_.size()),
    massTransferRates_(species_.size()),
    inertIndex_(species_[dict.lookup<word>("inertSpecie")]),
    active_(species_.size(), true)
{
    tmp<volScalarField> tYdefault;
    volScalarField YTot(volScalarField::New("YTot", mesh, 0.0));

    forAll(species_, i)
    {
        IOobject header
        (
            IOobject::groupName(species_[i], phaseName),
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
                        IOobject::groupName(species_[i], phaseName),
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
                word YdefaultName(IOobject::groupName("Ydefault", phaseName));

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
                        IOobject::groupName(species_[i], phaseName),
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
                    "massTransfer:" + IOobject::groupName(species_[i], phaseName),
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

Foam::multicomponentBlastThermo::~multicomponentBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multicomponentBlastThermo::initializeModels()
{
    if
    (
        mesh_.template foundObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", phaseName_)
        )
    )
    {
        alphaRhoPhiPtr_.set
        (
            &mesh_.template lookupObject<surfaceScalarField>
            (
                IOobject::groupName("alphaRhoPhi", phaseName_)
            )
        );
    }
    else if
    (
        mesh_.template foundObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", masterName_)
        )
    )
    {
        alphaRhoPhiPtr_.set
        (
            &mesh_.template lookupObject<surfaceScalarField>
            (
                IOobject::groupName("alphaRhoPhi", masterName_)
            )
        );
    }
    else
    {
        alphaRhoPhiPtr_.set
        (
            &mesh_.template lookupObject<surfaceScalarField>
            (
                IOobject::groupName("rhoPhi", masterName_)
            )
        );
    }

    if
    (
        mesh_.template foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", phaseName_)
        )
    )
    {
        alphaRhoPtr_.set
        (
            &mesh_.template lookupObject<volScalarField>
            (
                IOobject::groupName("alphaRho", phaseName_)
            )
        );
    }
    else if
    (
        mesh_.template foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", masterName_)
        )
    )
    {
        alphaRhoPtr_.set
        (
            &mesh_.template lookupObject<volScalarField>
            (
                IOobject::groupName("alphaRho", masterName_)
            )
        );
    }
    else
    {
        alphaRhoPtr_.set
        (
            &mesh_.template lookupObject<volScalarField>
            (
                IOobject::groupName("rho", masterName_)
            )
        );
    }
}


void Foam::multicomponentBlastThermo::solve()
{
    const dimensionedScalar& dT(mesh_.time().deltaT());
    volScalarField YTot(volScalarField::New("YTot", mesh_, 0.0));
    const surfaceScalarField& alphaRhoPhi(alphaRhoPhiPtr_());
    const volScalarField& alphaRho(alphaRhoPtr_());
    dimensionedScalar residualRho(dimDensity, 1e-10);

    forAll(Ys_, i)
    {
        if (active(i) && i != inertIndex_)
        {
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
                )/max(residualRho, alphaRho);
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


void Foam::multicomponentBlastThermo::addDelta
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
}

// ************************************************************************* //
