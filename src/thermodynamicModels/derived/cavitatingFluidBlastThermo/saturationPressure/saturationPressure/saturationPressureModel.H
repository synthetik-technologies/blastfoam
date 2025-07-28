/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::saturationPressureModel

Description
    Model to add additional energy based on saturationPressure reactions to increase the
    pressure of the JWL model.

SourceFiles
    saturationPressureModel.C
    newAfterburnModel.C

\*---------------------------------------------------------------------------*/

#ifndef saturationPressureModel_H
#define saturationPressureModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dimensionedScalar.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class saturationPressureModel Declaration
\*---------------------------------------------------------------------------*/

class saturationPressureModel
{
protected:

    //- Offset temperature
    dimensionedScalar Toffset_;


public:

    //- Runtime type information
    TypeName("saturationPressureModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            saturationPressureModel,
            dictionary,
            (const dictionary& dict),
            (dict)
        );

    // Constructor
    saturationPressureModel(const dictionary& dict);


    //- Destructor
    virtual ~saturationPressureModel();


    // Selectors

        static autoPtr<saturationPressureModel> New
        (
            const word& name,
            const dictionary& dict
        );


    // Member functions

        //- Return satuation pressure given a temperature
        virtual scalar pSat(const scalar T) const = 0;

        //- Return derivative of satuation pressure w.r.t temperature
        virtual scalar derivative(const scalar T) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
