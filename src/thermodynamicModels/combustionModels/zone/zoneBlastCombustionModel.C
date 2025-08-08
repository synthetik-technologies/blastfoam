/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "zoneBlastCombustionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(zone, 0);

    addToRunTimeSelectionTable(blastCombustionModel, zone, fluid);
    addNamedToRunTimeSelectionTable
    (
        blastCombustionModel,
        zone,
        fluid,
        zoneCombustion
    );

    addToRunTimeSelectionTable(blastCombustionModel, zone, solid);
    addNamedToRunTimeSelectionTable
    (
        blastCombustionModel,
        zone,
        solid,
        zoneCombustion
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::zone::filter
(
    const tmp<fvScalarMatrix>& tR
) const
{
    fvScalarMatrix& R = tR.ref();
    scalarField& Su = R.source();
    scalarField filteredField(Su.size(), 0);

    forAll(zoneNames_, zonei)
    {
        const labelList& cells = this->mesh().cellZones()[zoneNames_[zonei]];

        forAll(cells, i)
        {
            filteredField[cells[i]] = Su[cells[i]];
        }
    }

    Su = filteredField;

    if (R.hasDiag())
    {
        scalarField& Sp = R.diag();

        forAll(zoneNames_, zonei)
        {
            const labelList& cells =
                this->mesh().cellZones()[zoneNames_[zonei]];

            forAll(cells, i)
            {
                filteredField[cells[i]] = Sp[cells[i]];
            }
        }

        Sp = filteredField;
    }

    return tR;
}


template<class GeoField>
inline Foam::tmp<GeoField>
Foam::combustionModels::zone::filter
(
    const tmp<GeoField>& tS
) const
{
    scalarField& S = tS.ref();
    scalarField filteredField(S.size(), 0);

    forAll(zoneNames_, zonei)
    {
        const labelList& cells = this->mesh().cellZones()[zoneNames_[zonei]];

        forAll(cells, i)
        {
            filteredField[cells[i]] = S[cells[i]];
        }
    }

    S = filteredField;

    return tS;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::zone::zone
(
    const word& modelType,
    const multicomponentBlastThermo& thermo,
    const bool isFluid,
    const word& combustionProperties
)
:
    blastCombustionModel
    (
        modelType,
        thermo,
        combustionProperties
    ),
    combustionModelPtr_
    (
        isFluid
      ? blastCombustionModel::NewFluid
        (
            thermo,
            "zoneProperties"
        )
      : blastCombustionModel::NewSolid
        (
            thermo,
            "zoneProperties"
        )
    ),
    zoneNames_(this->coeffs().lookup("zones"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::zone::~zone()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::zone::correct()
{
    combustionModelPtr_->correct();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::combustionModels::zone::R(const label speciei) const
{
    return filter(combustionModelPtr_->R(speciei));
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::zone::R(volScalarField& Y) const
{
    return filter(combustionModelPtr_->R(Y));
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::zone::Qdot() const
{
    return filter(combustionModelPtr_->Qdot());
}


bool Foam::combustionModels::zone::read()
{
    if (blastCombustionModel::read())
    {
        combustionModelPtr_->read();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
