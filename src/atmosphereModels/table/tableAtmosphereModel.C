/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "tableAtmosphereModel.H"
#include "thermodynamicConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace atmosphereModels
{
    defineTypeNameAndDebug(table, 0);
    addToRunTimeSelectionTable(atmosphereModel, table, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmosphereModels::table::table
(
    const fvMesh& mesh,
    const dictionary& dict,
    const label zoneID
)
:
    atmosphereModel(mesh, dict, zoneID),
    pTable_(dict_.subDict("pTable"), "h", "p"),
    TTable_(dict_.subDict("TTable"), "h", "T"),
    correct_(dict_.lookupOrDefault("correct", false))
{
    const_cast<dictionary&>(dict_).set
    (
        "pRef",
        pTable_.lookup(gMin(h_))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::atmosphereModels::table::~table()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmosphereModels::table::createAtmosphere
(
    fluidBlastThermo& thermo
) const
{
    volScalarField& p = thermo.p();
    volScalarField& T = thermo.T();

    // Optional setting of only some cells
    labelList cells
    (
        zoneID_ >= 0
      ? labelList(mesh_.cellZones()[zoneID_])
      : identity(mesh_.nCells())
    );

    // Create a hash set for easier searching of selected cells
    labelHashSet cSet(cells);

    // Set the internal values
    forAll(cells, i)
    {
        p[cells[i]] = pTable_.lookup(h_[cells[i]]);
        T[cells[i]] = TTable_.lookup(h_[cells[i]]);
    }

    // The the boundary values that have their owner face included in the set
    volScalarField::Boundary& bp = p.boundaryFieldRef();
    volScalarField::Boundary& bT = T.boundaryFieldRef();
    forAll(bp, patchi)
    {
        const labelList& fCells = bp[patchi].patch().faceCells();
        forAll(fCells, facei)
        {
            if (cSet.found(fCells[facei]))
            {
                bp[patchi][facei] =
                    pTable_.lookup(h_.boundaryField()[patchi][facei]);
                bT[patchi][facei] =
                    TTable_.lookup(h_.boundaryField()[patchi][facei]);
            }
        }
    }

    // Correct boundary conditions
    p.correctBoundaryConditions();
    T.correctBoundaryConditions();

    // Correct density
    thermo.updateRho(p);

    // Correct of thermodynamic variables
    thermo.correct();

    // Equalibriate the pressure field
    if (correct_)
    {
        hydrostaticInitialisation(thermo);
    }
}

// ************************************************************************* //
