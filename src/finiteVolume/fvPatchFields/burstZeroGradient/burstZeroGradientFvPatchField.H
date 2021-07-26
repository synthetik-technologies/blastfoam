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
    Foam::burstZeroGradientFvPatchVectorField

Description
    This velocity boundary condition simulates the opening of a ZeroGradient due
    to local flow conditions, by merging the behaviours of wall and cyclic
    conditions.  The ZeroGradient joins two mesh regions, where the open fraction
    determines the interpolation weights applied to each cyclic- and
    neighbour-patch contribution.

    We determine whether the ZeroGradient is opening or closing from the sign of
    the net force across the ZeroGradient, from which the ZeroGradient open fraction is
    updated using:

        \f[
            x = x_{old} + sign(F_{net})\frac{dt}{DT}
        \f]

    where

    \vartable
        x       | ZeroGradient open fraction [0-1]
        x_{old} | ZeroGradient open fraction on previous evaluation
        dt      | simulation time step
        DT      | time taken to open the ZeroGradient
        F_{net} | net force across the ZeroGradient
    \endvartable

    The open fraction is then applied to scale the patch areas.

Usage
    \table
        Property     | Description             | Required    | Default value
        p            | pressure field name     | no          | p
        cyclicPatch  | cyclic patch name      | yes         |
        orientation  | 1 or -1 used to switch flow direction | yes|
        openFraction | current opatch open fraction [0-1]| yes |
        openingTime  | time taken to open the ZeroGradient | yes |
        maxOpenFractionDelta | max open fraction change per timestep | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            burstZeroGradient;
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
    burstZeroGradientFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef burstZeroGradientFvPatchVectorField_H
#define burstZeroGradientFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class burstZeroGradientFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class burstZeroGradientFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Name of the pressure field used to calculate the force
        //  on the active ZeroGradient
        word pName_;

        //- Name of the cyclic patch used when the active ZeroGradient is open
        word cyclicPatchName_;

        //- Index of the cyclic patch used when the active ZeroGradient is open
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
    TypeName("burstZeroGradient");


    // Constructors

        //- Construct from patch and internal field
        burstZeroGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        burstZeroGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given burstZeroGradientFvPatchField
        //  onto a new patch
        burstZeroGradientFvPatchField
        (
            const burstZeroGradientFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        burstZeroGradientFvPatchField
        (
            const burstZeroGradientFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new burstZeroGradientFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        burstZeroGradientFvPatchField
        (
            const burstZeroGradientFvPatchField<Type>&,
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
                new burstZeroGradientFvPatchField<Type>(*this, iF)
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
#include "burstZeroGradientFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
