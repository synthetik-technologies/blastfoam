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
    Foam::globalTemperatureCoupledFvPatchScalarField

Description
    Mixed boundary condition for temperature and radiation heat transfer
    to be used for in multiregion cases. Optional thin thermal layer
    resistances can be specified through thicknessLayers and kappaLayers
    entries.
    The thermal conductivity \c kappa can either be retrieved from various
    possible sources.

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | T
        qrNbr      | name of the radiative flux in the nbr region | no | none
        qr         | name of the radiative flux in this region | no | none
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivites per layer [W/m/K] | no |
    \endtable
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            globalTemperatureCoupled;
        Tnbr            T;
        qrNbr           qr; // or none. Name of qr field on neighbour region
        qr              qr; // or none. Name of qr field on local region
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);
        value           uniform 300;
    }
    \endverbatim

See also

SourceFiles
    globalTemperatureCoupledFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef globalTemperatureCoupledFvPatchScalarField_H
#define globalTemperatureCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "globalPolyBoundaryMesh.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class globalTemperatureCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class globalTemperatureCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Reference to the global boundaryMesh
        const globalPolyBoundaryMesh& globalBoundary_;

        //- Name of field on the neighbour region
        const word TnbrName_;

         //- Name of the radiative heat flux in the neighbor region
        const word qrNbrName_;

        //- Name of the radiative heat flux in local region
        const word qrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        //- Total contact resistance
        scalar contactRes_;

        //- Unmapped temperature
        scalar unmappedT_;


    // Private member functions

        //- Given patch temperature calculate corresponding K field
        static tmp<scalarField> kappa(const fvPatchScalarField& Tp);

public:

    //- Runtime type information
    TypeName("globalTemperatureCoupled");


    // Constructors

        //- Construct from patch and internal field
        globalTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        globalTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  globalTemperatureCoupledFvPatchScalarField onto a
        //  new patch
        globalTemperatureCoupledFvPatchScalarField
        (
            const globalTemperatureCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        globalTemperatureCoupledFvPatchScalarField
        (
            const globalTemperatureCoupledFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new globalTemperatureCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
