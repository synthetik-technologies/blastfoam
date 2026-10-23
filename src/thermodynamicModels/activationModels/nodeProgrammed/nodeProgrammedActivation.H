/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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
    Foam::activationModels::nodeProgrammedActivation

Description
    Node/vertex-based programmed burn activation model.
    Computes cell arrival time extrema (tStart and tEnd) from cell vertices
    and smoothly ramps the burn progress variable lambda across the cell:
    lambda = (t - tStart) / (tEnd - tStart).

SourceFiles
    nodeProgrammedActivation.C

\*---------------------------------------------------------------------------*/

#ifndef nodeProgrammedActivation_H
#define nodeProgrammedActivation_H

#include "activationModel.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace activationModels
{

/*---------------------------------------------------------------------------*\
                   Class nodeProgrammedActivation Declaration
\*---------------------------------------------------------------------------*/

class nodeProgrammedActivation
:
    public activationModel
{
public:
    //- Burn models
        enum burnModel
        {
            PROGRAMMED,     // Node-based geometric burn
            PROGRAMMEDBETA  // Node-based programmed + compression burn
        };

        static const NamedEnum<burnModel, 2> burnModelNames_;

private:
// Private data

        //- Const reference to density
        const volScalarField& rho_;

        //- Name of pressure field
        word pName_;

        //- Const reference to pressure
        const volScalarField& p_;

        //- Minimum pressure for compression ignition
        dimensionedScalar pMin_;

        //- Pointer to volume fraction field
        const volScalarField* alphaPtr_;

        //- Residual volume fraction
        scalar residualAlpha_;

        //- Detonation velocity
        dimensionedScalar vDet_;

        //- Chapman-Jouguet pressure
        dimensionedScalar Pcj_;

        //- Reference density
        dimensionedScalar rho0_;

        //- Chapman-Jouguet specific volume
        dimensionedScalar Vcj_;

        //- Burn model
        burnModel model_;

        //- Cell entry time (minimum vertex arrival time)
        volScalarField::Internal tStart_;

        //- Cell exit time (maximum vertex arrival time)
        volScalarField::Internal tEnd_;

        //- Add the delay of the detonation points to the arrival times
        Switch useDelay_;

        //- Time index of the last calculation of the arrival times
        label timeIndex_;

        //- Calculate the cell arrival times from the current mesh points
        //  (on construction and after topology changes, e.g. refinement,
        //  since arrival times mapped between cell levels are wrong)
        void calcArrivalTimes();

        //- Return the time rate of change of lambda
        virtual tmp<volScalarField> delta() const;

        //- Correct lambda field after advection
        virtual void correct(volScalarField& lambda);


public:

    //- Runtime type information
    TypeName("nodeProgrammed");

    // Constructor
    nodeProgrammedActivation
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );

    //- Destructor
    virtual ~nodeProgrammedActivation();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
