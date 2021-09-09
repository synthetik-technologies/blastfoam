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
    Foam::FieldSetTypes::MassIntegrated

Description
    Given a mass integrated quantity, compute a field value

Usage:
    volScalarFieldMassIntegrated alpha.particles 10.23


SourceFiles
    MassIntegratedFieldSetType.C
    MassIntegratedFieldSetType.C

\*---------------------------------------------------------------------------*/

#ifndef MassIntegratedFieldSetType_H
#define MassIntegratedFieldSetType_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldSetType.H"
#include "blastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldSetTypes
{

/*---------------------------------------------------------------------------*\
                          Class MassIntegrated Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class MassIntegrated
:
    public FieldSetType<Type, Patch, Mesh>
{
    //- Name of the thermo phase
    word phaseName_;

    //- Value to set
    Type value_;

    //- Dummy field for solid phase pressures
    autoPtr<volScalarField> psPtr_;

    //- Thermodynamic model
    const blastThermo& thermo_;

    //- Set the selected cells
    virtual void setField();

    //- Read the phase name of the field
    word readPhaseName(Istream& is, const word& fieldName) const;

    //- Return the correct dictionary to create the thermo model from
    dictionary thermoDict(const fvMesh& mesh, const word& phaseName) const;

    //- Either create of lookup thermo so it only needs to be created once
    const blastThermo& lookupOrConstructThermo
    (
        const fvMesh& mesh,
        const word& phaseName
    ) const;


public:

    TypeName("MassIntegrated");

    // Constructor
    MassIntegrated
    (
        const fvMesh& mesh,
        const word& fieldName,
        const labelList& selectedCells,
        Istream& is,
        const bool write
    );

    //- Destructor
    virtual ~MassIntegrated();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FieldSetTypes
} // End namespace Foam

#ifdef NoRepository
    #include "MassIntegratedFieldSetType.C"
#endif

#endif
// ************************************************************************* //
