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
    Foam::interfacialPressureModels::massAveraged

Description
    One phase is used to determine the interfacial pressure

SourceFiles
    massAveragedPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef massAveragedInterfacialPressure_H
#define massAveragedInterfacialPressure_H

#include "interfacialPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfacialPressureModels
{

/*---------------------------------------------------------------------------*\
                          Class massAveraged Declaration
\*---------------------------------------------------------------------------*/

class massAveraged
:
    public interfacialPressureModel
{
    // Private data

        //- List of phase names used to determine interfacial pressure
        wordList phaseNames_;


public:

    //- Runtime type information
    TypeName("massAveraged");


    // Constructors

        //- Construct from components
        massAveraged
        (
            const dictionary& dict,
            const phaseModelList& phaseModels
        );


    //- Destructor
    virtual ~massAveraged();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        virtual tmp<volScalarField> Pi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfacialPressureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
