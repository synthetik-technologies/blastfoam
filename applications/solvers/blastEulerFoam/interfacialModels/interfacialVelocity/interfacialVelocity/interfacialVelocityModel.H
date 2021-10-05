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
    Foam::interfacialVelocityModel

Description

SourceFiles
    interfacialVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfacialVelocityModel_H
#define interfacialVelocityModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseModel.H"
#include "phasePair.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class interfacialVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class interfacialVelocityModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Reference to list of phaseModels
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("interfacialVelocityModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            interfacialVelocityModel,
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
        interfacialVelocityModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfacialVelocityModel();


    // Selectors

        static autoPtr<interfacialVelocityModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Return the interfacial Velocity
        virtual tmp<volVectorField> UI() const = 0;

        //- Return the interfacial volumetric flux
        virtual tmp<surfaceScalarField> phiI() const = 0;

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
