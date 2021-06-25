/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-06-2021 Jeff Heylmun     | Created class for solid model
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
    Foam::linearElastic

Description
    Linear elastic solid model

SourceFiles
    linearElastic.C
    newlinearElastic.C

\*---------------------------------------------------------------------------*/

#ifndef linearElastic_H
#define linearElastic_H

#include "solidModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "solidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                           Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class linearElastic
:
    public solidModel
{
    // Private data

        //- Mechanical properties dictionary
        IOdictionary mechanicalProperties_;

        //- Is plane stress used
        Switch planeStress_;

        volScalarField E_;

        volScalarField nu_;

        volScalarField mu_;

        volScalarField lambda_;

        //- Is thermal stress used
        Switch thermalStress_;

        autoPtr<volScalarField> alphaPtr_;

        //- Compact normal stress
        Switch compactNormalStress_;

        //- Stress field
        //  This is the engineering stress for linear geometry approaches, and
        //  true (Cauchy) stress for nonlinear geometry (large strain)
        //  approaches
        volSymmTensorField sigmaD_;

        volVectorField divSigmaExp_;

        volScalarField sigmaEq_;




public:

    //- Runtime type information
    TypeName("linearElastic");

    // Constructors

        //- Construct from components
        linearElastic(fvMesh& mesh);

    // Destructor

        virtual ~linearElastic();


    // Member Functions


        // Edit

            //- Evolve the solid model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
