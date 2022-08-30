/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::thermoBasePatchField

Description
    Basic class to return thermodynamic properties

SourceFiles
    thermoBasePatchField.C

See also
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef thermoBasePatchField_H
#define thermoBasePatchField_H

#include "fixedValueFvPatchFields.H"
#include "UautoPtr.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class thermoBasePatchField Declaration
\*---------------------------------------------------------------------------*/

class thermoBasePatchField
{
protected:

    // Private Data

        //- Reference to the patch
        const fvPatch& patch_;

        //- Name of the velocity field
        word phaseName_;

        //- Pointer to the thermo class
        mutable UautoPtr<const fluidThermo> thermoPtr_;


public:

    //- Runtime type information
    TypeName("thermoBase");


    // Constructors

        //- Construct from patch
        thermoBasePatchField(const fvPatch&);

        //- Construct from patch and internal field
        thermoBasePatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch and dictionary
        thermoBasePatchField
        (
            const fvPatch&,
            const dictionary&
        );

        //- Construct from patch, internal field, and dictionary
        thermoBasePatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );


    // Destructor

        virtual ~thermoBasePatchField();


    // Member Functions

        // Access

            //- Return the fluidThermo of this phase
            const fluidThermo& thermo() const;

            //- Return the psi field
            tmp<scalarField> psi() const;

            //- Return the gamma field
            tmp<scalarField> gamma() const;

            //- Return the gamma field given a temperature
            tmp<scalarField> gamma(const scalarField& T) const;

            //- Return the speed of sound
            tmp<scalarField> speedOfSound() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
