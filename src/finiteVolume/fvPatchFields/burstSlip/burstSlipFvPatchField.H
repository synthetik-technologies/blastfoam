/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::burstSlipFvPatchVectorField

Description
    This velocity boundary condition simulates the opening of a Slip due
    to local flow conditions, by merging the behaviours of wall and cyclic
    conditions.  The Slip joins two mesh regions, where the open fraction
    determines the interpolation weights applied to each cyclic- and
    neighbour-patch contribution.

    We determine whether the Slip is opening or closing from the sign of
    the net force across the Slip, from which the Slip open fraction is
    updated using:

        \f[
            x = x_{old} + sign(F_{net})\frac{dt}{DT}
        \f]

    where

    \vartable
        x       | Slip open fraction [0-1]
        x_{old} | Slip open fraction on previous evaluation
        dt      | simulation time step
        DT      | time taken to open the Slip
        F_{net} | net force across the Slip
    \endvartable

    The open fraction is then applied to scale the patch areas.

Usage
    \table
        Property     | Description             | Required    | Default value
        p            | pressure field name     | no          | p
        cyclicPatch  | cyclic patch name      | yes         |
        orientation  | 1 or -1 used to switch flow direction | yes|
        openFraction | current opatch open fraction [0-1]| yes |
        openingTime  | time taken to open the Slip | yes |
        maxOpenFractionDelta | max open fraction change per timestep | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            burstSlip;
        p               p;
        cyclicPatch     cyclic1;
        orientation     1;
        openFraction    0.2;
        openingTime     5.0;
        maxOpenFractionDelta 0.1;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::cyclicFvPatchField

SourceFiles
    burstSlipFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef burstSlipFvPatchSlipField_H
#define burstSlipFvPatchSlipField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class burstSlipFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class burstSlipFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Name of the pressure field used to calculate the force
        //  on the active Slip
        word pName_;

        //- Name of the cyclic patch used when the active Slip is open
        word cyclicPatchName_;

        //- Index of the cyclic patch used when the active Slip is open
        label cyclicPatchLabel_;

        //- Burst pressure
        scalar pBurst_;

        //- Burst impulse
        scalar burstImpulse_;

        //- Is impulse used
        Switch useImpulse_;

        //- Impulse field
        scalarField impulse_;

        //- Do individual faces burst
        Switch partialBurst_;

        //- Intact field
        scalarField intact_;

        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("burstSlip");


    // Constructors

        //- Construct from patch and internal field
        burstSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        burstSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given burstSlipFvPatchField
        //  onto a new patch
        burstSlipFvPatchField
        (
            const burstSlipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        burstSlipFvPatchField
        (
            const burstSlipFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        burstSlipFvPatchField
        (
            const burstSlipFvPatchField<Type>&,
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
                new burstSlipFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchField<Type>&, const labelList&);


    // Evaluation functions

        //- Return neighbour coupled internal cell data
        tmp<Field<Type>> patchNeighbourField() const;

        //- Return neighbour coupled internal cell data
        tmp<Field<Type>> truePatchNeighbourField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "burstSlipFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
