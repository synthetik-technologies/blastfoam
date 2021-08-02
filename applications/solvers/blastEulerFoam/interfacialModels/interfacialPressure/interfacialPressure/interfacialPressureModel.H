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
    Foam::interfacialPressureModel

Description

SourceFiles
    interfacialPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfacialPressureModel_H
#define interfacialPressureModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseModel.H"
#include "runTimeSelectionTables.H"
#include "phasePair.H"
#include "fluidBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class interfacialPressureModel Declaration
\*---------------------------------------------------------------------------*/

class interfacialPressureModel
:
    public regIOobject
{
protected:

    //- Const reference to phase pair
    const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("interfacialPressureModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            interfacialPressureModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static data members


    // Constructors

        //- Construct froma dictionary and a phase pair
        interfacialPressureModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfacialPressureModel();


    // Selectors

        static autoPtr<interfacialPressureModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Return the interfacial pressure
        virtual tmp<volScalarField> Pi() const = 0;

        //- Return the interfacial pressure for a cell
        virtual scalar Pi(const label celli) const = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const
        {
            return os.good();
        }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
