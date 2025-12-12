/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
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
    Foam::dragModels::dispersedHeatTransferModel

Description

SourceFiles
    dispersedHeatTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedHeatTransferModel_H
#define dispersedHeatTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "heatTransferModel.H"

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                Class dragModels::dispersedHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedHeatTransferModel
:
    public heatTransferModel
{
protected:

    // Protected data

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("dispersedHeatTransferModel");


    // Constructors

        // Construct with residual constants
        dispersedHeatTransferModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~dispersedHeatTransferModel();


    // Member Functions

        //- Drag coefficient

        //- Return the phase-intensive heat transfer coefficient Ki
        //  used in the enthalpy equations
        //    ddt(alpha1*rho1*ha) + ... = ... alphad*K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... alphad*K*(Tb - Ta)
        virtual tmp<volScalarField> Ki() const = 0;

        //- Return the phase-intensive heat transfer coefficient Ki
        //  used in the enthalpy equations
        //    ddt(alpha1*rho1*ha) + ... = ... alphad*K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... alphad*K*(Tb - Ta)
        virtual scalar cellKi(const label) const = 0;

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        //  with a specified residual volume fraction
        virtual tmp<volScalarField> stabK() const;

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        //  with a specified residual volume fraction
        virtual scalar cellStabK(const label celli) const;

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        //  with a specified residual volume fraction
        virtual tmp<volScalarField> K() const;

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        //  with a specified residual volume fraction
        virtual scalar cellK(const label celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
