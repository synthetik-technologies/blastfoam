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

Class
    Foam::pressureWaveTransmissiveFvPatchField

Description
    This boundary condition provides a wave transmissive outflow condition,
    based on solving DDt(W, field) = 0 at the boundary \c W is the wave velocity
    and speed of sound field is the field to which this boundary condition is
    applied.

    The wave speed is calculated using:

        \f[
            w_p = \frac{\phi_p}{|Sf|} + c
        \f]

    where

    \vartable
        w_p     | patch wave speed
        \phi_p  | patch face flux
        Sf      | patch face area vector
        c       | speed of sound
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        phaseName    | phase name              | no          | word::null
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            waveTransmissive;
        phi             phi;
    }
    \endverbatim

See also
    Foam::advectiveFvPatchField

SourceFiles
    pressureWaveTransmissiveFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureWaveTransmissiveFvPatchField_H
#define pressureWaveTransmissiveFvPatchField_H

#include "advectiveFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class pressureWaveTransmissiveFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pressureWaveTransmissiveFvPatchField
:
    public advectiveFvPatchField<Type>
{

    // Private Data

        //- Name of the phase
        word phaseName_;


public:

    //- Runtime type information
    TypeName("pressureWaveTransmissive");


    // Constructors

        //- Construct from patch and internal field
        pressureWaveTransmissiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureWaveTransmissiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureWaveTransmissiveFvPatchField
        //  onto a new patch
        pressureWaveTransmissiveFvPatchField
        (
            const pressureWaveTransmissiveFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        pressureWaveTransmissiveFvPatchField
        (
            const pressureWaveTransmissiveFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new pressureWaveTransmissiveFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        pressureWaveTransmissiveFvPatchField
        (
            const pressureWaveTransmissiveFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new pressureWaveTransmissiveFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions


        // Evaluation functions

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pressureWaveTransmissiveFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
