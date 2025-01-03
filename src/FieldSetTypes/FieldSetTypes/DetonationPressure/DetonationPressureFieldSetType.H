/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

Class
    Foam::FieldSetTypes::DetonationPressure

Description
    Given phase name, detonation energy, and reference density, set the pressure

Usage:
    volScalarFieldDetonationPressure p tnt 7e9 1630


SourceFiles
    DetonationPressureFieldSetType.C
    DetonationPressureFieldSetType.C

\*---------------------------------------------------------------------------*/

#ifndef DetonationPressureFieldSetType_H
#define DetonationPressureFieldSetType_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldSetType.H"
#include "fluidBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldSetTypes
{

/*---------------------------------------------------------------------------*\
                          Class DetonationPressure Declaration
\*---------------------------------------------------------------------------*/

class DetonationPressure
:
    public VolFieldSetType<scalar>
{
    //- Name of the thermo phase
    word phaseName_;

    //- Specific detonation energy
    scalar e0_;

    //- Calculated detonation pressure
    volScalarField* pDetPtr_;


    //- Read the phase name of the field
    word readPhaseName(Istream& is, const word& fieldName) const;

    //- Return the correct dictionary to create the thermo model from
    dictionary thermoDict(const fvMesh& mesh, const word& phaseName) const;

    //- Either create of lookup thermo so it only needs to be created once
    fluidBlastThermo& lookupOrConstructThermo
    (
        const fvMesh& mesh,
        const word& phaseName
    ) const;

    //- Return the subset of values given indices and positions
    virtual void setGeoField
    (
        const labelList& indices,
        const UIndirectList<vector>& pts,
        UIndirectList<scalar>& f,
        const label patchi
    );


public:

    TypeName("DetonationPressure");

    // Constructor
    DetonationPressure
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& fieldName,
        const labelList& selectedCells,
        Istream& is,
        const bool write
    );

    //- Destructor
    virtual ~DetonationPressure();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FieldSetTypes
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
