/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
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

#include "MixtureThermoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo>
Foam::MixtureThermoModel<BasicThermo>::MixtureThermoModel
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
    }
    forAll(Ys_, i)
    {
        Ys_[i] /= YTot;
        Ys_[i].correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo>
Foam::MixtureThermoModel<BasicThermo>::~MixtureThermoModel()
{}

// ************************************************************************* //
