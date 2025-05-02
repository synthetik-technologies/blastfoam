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
    functionObjects::hotCylinderAnalyticalSolution

Description
    FunctionObject that generates the analytical temperature (analyticalT) and
    stress fields (analyticalRadialStress, analyticalHoopStress) for the case of
    a thermally-stressed pipe/cylinder.
    The solution is taken from Timoshenko, Theory of Elasticity.

    The analytical solution requires the following inputs:
        - rInner : inner pipe radius
        - rOuter : outer pipe radius
        - TInner : temperature on the inner pipe surface
        - TOuter : temperature on the outer pipe surface
        - E : Young's modulus
        - nu : Poisson's ratio
        - alpha : coefficient of linear thermal expansion

    It is assumed that the inner and outer pipe surfaces are traction-free.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    hotCylinderAnalyticalSolution.C

\*---------------------------------------------------------------------------*/

#ifndef hotCylinderAnalyticalSolution_H
#define hotCylinderAnalyticalSolution_H

#include "fvMeshFunctionObject.H"
#include "dictionary.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class hotCylinderAnalyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class hotCylinderAnalyticalSolution
:
    public fvMeshFunctionObject
{
    // Private data

        //- Inner pipe radius
        dimensionedScalar rInner_;

        //- Outer pipe radius
        dimensionedScalar rOuter_;

        //- Inner pipe temperature
        dimensionedScalar TInner_;

        //- Outer pipe temperature
        dimensionedScalar TOuter_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Coefficient of linear thermal expansion
        dimensionedScalar alpha_;


public:

    //- Runtime type information
    TypeName("hotCylinderAnalyticalSolution");


    // Constructors

        //- Construct from components
        hotCylinderAnalyticalSolution
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        virtual wordList fields() const
        {
            return wordList::null();
        };

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
