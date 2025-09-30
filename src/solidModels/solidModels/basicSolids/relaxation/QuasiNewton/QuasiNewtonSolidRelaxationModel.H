/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    QuasiNewtonSolidRelaxationModel

Description
    Quasi-Newton solid relaxation

SourceFiles
    QuasiNewtonSolidRelaxationModel.C

\*---------------------------------------------------------------------------*/

#ifndef QuasiNewtonSolidRelaxationModel_H
#define QuasiNewtonSolidRelaxationModel_H

#include "solidRelaxationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class QuasiNewtonSolidRelaxationModel Declaration
\*---------------------------------------------------------------------------*/

class QuasiNewtonSolidRelaxationModel
:
    public solidRelaxationModel
{
    // Private data

        // Quasi-Newton fields

            //- Restart frequency for the Quasi-Newton method
            label restartFreq_;

            //- Quasi-Newton method input vector field
            DynamicList<vectorField> V_;

            //- Quasi-Newton method output vector field
            DynamicList<vectorField> W_;

            //- Quasi-Newton method time indices for the input-ouput fields
            DynamicList<scalar> T_;

            //- Reference D field used for Quasi-Newton method
            volVectorField DRef_;

            //- Reference D (without relaxation) field used for Quasi-Newton method
            volVectorField unrelaxedDRef_;

    // Protected member functions

        //- Read parameters
        virtual void readDict();


public:

    // Constructors

        //- Construct from components
        QuasiNewtonSolidRelaxationModel(solidModel& solid);


    // Destructor
    virtual ~QuasiNewtonSolidRelaxationModel();


    //- Apply fixed or Aitken's adaptive under-relaxation to the field
    void relaxField(volVectorField& D, const label iCorr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
