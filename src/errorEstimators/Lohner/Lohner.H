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
    Foam::errorEstimators::Lohner

Description
    Estimates error based on second derivative of a specified field

    References:
    \verbatim
        R. Löhner (1987).
        An adaptive finite element scheme for transient problems in
 CFD
        Computer Methods in Applied Mechanics and Engineering, 61, 323-338.
    \endverbatim

SourceFiles
    Lohner.C

\*---------------------------------------------------------------------------*/

#ifndef Lohner_H
#define Lohner_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "errorEstimator.H"

namespace Foam
{
namespace errorEstimators
{

/*---------------------------------------------------------------------------*\
                           Class Lohner Declaration
\*---------------------------------------------------------------------------*/

class Lohner
:
    public errorEstimator
{
// Private data

        //- Field name
        const word fieldName_;

        //- Coefficient
        const scalar epsilon_;


public:

    //- Runtime type information
    ClassName("Lohner");

    // Constructor
    Lohner(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~Lohner();


    // Member Functions

        //- Gamma function
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorEstimators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
