/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "basicSpecieBlastMixture.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicSpecieBlastMixture, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::basicSpecieBlastMixture::normaliseMassFractions()
{
    if (!species_.size())
    {
        return;
    }

    tmp<volScalarField> tYt
    (
        volScalarField::New
        (
            IOobject::groupName("Yt", Y_[0].group()),
            Y_[0].mesh(),
            0.0,
            calculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& Yt = tYt.ref();

    forAll(Y_, i)
    {
        Yt += Y_[i];
    }

    if (mag(min(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << species() << nl
            << "Min(Sum of mass fractions):
            << "    InternalField = " << gMin(Yt()) << nl;
            << "    BoundaryField = " << gMin(Yt.boundaryField()) << endl
            << exit(FatalError);
    }

    forAll(Y_, i)
    {
        Y_[i] /= Yt;
    }
}


void Foam::basicSpecieBlastMixture::correctMassFractions()
{
    if (!species_.size())
    {
        return;
    }

    if (defaultSpeciei_ < 0)
    {
        tmp<volScalarField> tYt
        (
            volScalarField::New
            (
                IOobject::groupName("Yt", Y_[0].group()),
                Y_[0].mesh(),
                0.0,
                calculatedFvPatchScalarField::typeName
            )
        );
        volScalarField& Yt = tYt.ref();

        forAll(Y_, i)
        {
            Yt += Y_[i];
        }

        if (mag(min(Yt).value()) < rootVSmall)
        {
            FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << species() << nl
            << "Min(Sum of mass fractions):
            << "    InternalField = " << gMin(Yt()) << nl;
            << "    BoundaryField = " << gMin(Yt.boundaryField()) << endl
            << exit(FatalError);
        }
        forAll(Y_, i)
        {
            Y_[i] /= Yt;
        }
    }
    else
    {
        volScalarField Ysum
        (
            volScalarField::New
            (
                IOobject::groupName("Ysum", Y_[0].group()),
                mesh(),
                0.0,
                calculatedFvPatchScalarField::typeName
            )
        );

        forAll(Y_, i)
        {
            if (i != defaultSpeciei_)
            {
                Ysum += Y_[i];
            }
        }

        Y()[defaultSpeciei_] = 1.0 - Ysum;
        Y()[defaultSpeciei_].max(0);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSpecieBlastMixture::basicSpecieBlastMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    phaseName_(phaseName),
    species_(specieNames),
    defaultSpecieName_
    (
        (species_.size() && phaseName.empty())
      ? thermoDict.lookupBackwardsCompatible<word>
        (
            {"defaultSpecie", "inertSpecie"}
        )
      : word::null
    ),
    defaultSpeciei_
    (
        species_.size()
     && species_.found(defaultSpecieName_)
      ? species_[defaultSpecieName_]
      : -1
    ),
    active_(species_.size(), true),
    Y_(species_.size())
{
    if (species_.size() && phaseName.empty() && defaultSpeciei_ == -1)
    {
        FatalIOErrorInFunction(thermoDict)
            << "default specie " << defaultSpecieName_
            << " not found in available species " << species_
            << exit(FatalIOError);
    }

    // Read the species' mass fractions
    tmp<volScalarField> tYdefault;
    forAll(species_, i)
    {
        typeIOobject<volScalarField> header
        (
            IOobject::groupName(species_[i], phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ
        );

        if (header.headerOk())
        {
            // Read the mass fraction field
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], phaseName),
                        mesh.time().name(),
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
                const word YdefaultName
                (
                    IOobject::groupName("Ydefault", phaseName)
                );

                typeIOobject<volScalarField> timeIO
                (
                    YdefaultName,
                    mesh.time().name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                typeIOobject<volScalarField> constantIO
                (
                    YdefaultName,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                typeIOobject<volScalarField> time0IO
                (
                    YdefaultName,
                    Time::timeName(0),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (timeIO.headerOk())
                {
                    tYdefault = new volScalarField(timeIO, mesh);
                }
                else if (constantIO.headerOk())
                {
                    tYdefault = new volScalarField(constantIO, mesh);
                }
                else
                {
                    tYdefault = new volScalarField(time0IO, mesh);
                }
            }

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tYdefault()
                )
            );
        }

        volScalarField& Y = Y_[i];
        volScalarField::Boundary& bY = Y.boundaryFieldRefNoStoreOldTimes();
        forAll(bY, patchi)
        {
            if (!bY[patchi].fixesValue())
            {
                bY[patchi] = bY[patchi].patchInternalField();
            }
        }
    }

    normaliseMassFractions();
}



// ************************************************************************* //
