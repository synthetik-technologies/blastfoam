/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::surfaceTemperatureModels::isoThermal

Description
    Surface temperature is equal to the cell centered temperature

SourceFiles
    isoThermalSurfaceTemperatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef isoThermalSurfaceTemperatureModel_H
#define isoThermalSurfaceTemperatureModel_H

#include "surfaceTemperatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTemperatureModels
{

/*---------------------------------------------------------------------------*\
                          Class isoThermal Declaration
\*---------------------------------------------------------------------------*/

class isoThermal
:
    public surfaceTemperatureModel
{
public:

    //- Runtime type information
    TypeName("isoThermal");


    // Constructors

        //- Construct from components
        isoThermal
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~isoThermal();


    // Member Functions

        //- Return the surface temperature
        virtual tmp<volScalarField> Ts() const;
        
        //- Update before each ode step
        virtual void update()
        {};

        //- Update after ode integration
        virtual void postUpdate()
        {};

        //- Solve sub-step stepi
        virtual void solve()
        {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTemperatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
