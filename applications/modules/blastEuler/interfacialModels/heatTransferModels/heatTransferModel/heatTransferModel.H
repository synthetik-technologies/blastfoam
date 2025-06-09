/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun:    Added cell based returns
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::heatTransferModel

Description

SourceFiles
    heatTransferModel.C
    heatTransferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferModel_H
#define heatTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "BlendedInterfacialModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                        Class heatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferModel
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("heatTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        //- Construct froma dictionary and a phase pair
        heatTransferModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~heatTransferModel();


    // Selectors

        static autoPtr<heatTransferModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        //  with a specified residual volume fraction
        virtual tmp<volScalarField> K() const = 0;

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        //  with a specified residual volume fraction
        virtual scalar cellK(const label celli) const = 0;
};


/*---------------------------------------------------------------------------*\
                      Class blendedDragModel Declaration
\*---------------------------------------------------------------------------*/

class blendedHeatTransferModel
:
    public BlendedInterfacialModel<heatTransferModel>
{
public:

    typedef heatTransferModel modelType;

    // Constructors

        //- Inherit base class constructors
        using BlendedInterfacialModel<heatTransferModel>::
            BlendedInterfacialModel;


    // Member Functions

        //- Return the drag coefficient K
        tmp<volScalarField> K() const;

        //- Return the blended cell force coefficient
        scalar cellK(const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
