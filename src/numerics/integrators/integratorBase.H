/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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
    Foam::Integrator

Description
    Generic class for integrating a given univariate equation.

SourceFiles
    integratorBase.C

\*---------------------------------------------------------------------------*/

#ifndef integratorBase_H
#define integratorBase_H

#include "typeInfo.H"
#include "adaptiveTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class integratorBase Declaration
\*---------------------------------------------------------------------------*/

class integratorBase
{
protected:
    // Protected data

        //- Is the integration adaptive
        bool adaptive_;

        //- Tracked number of function evaluations
        mutable label evals_;


public:

    //- Runtime type information
    TypeName("integratorBase");

    // Constructors

        //- Construct for from dictinoary
        integratorBase(const dictionary& dict);


    //- Destructor
    virtual ~integratorBase()
    {}


    // Member Functions

        //- Return the number of function evaluations
        label nEvals() const
        {
            return evals_;
        }

        //- Access the number of steps
        virtual void setAdaptive()
        {
            adaptive_ = true;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
