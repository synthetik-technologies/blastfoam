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
    Foam::timeIntegrators::RK3SSP

Description
    Second order, strong stability preserving Runge-Kutta method

    References:
    \verbatim
        Spiteri, R.J., Ruuth, S.J. (2003).
        A New Class of Optimal Higher-order Strong-Stability-Preserving
        Time Discretization Methods
        SIAM Journal of Numerical Analysis, 40(2), 469-491.
    \endverbatim

SourceFiles
    RK3SSPTimeIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef RK3SSPTimeIntegrator_H
#define RK3SSPTimeIntegrator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrator.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class RK3SSP Declaration
\*---------------------------------------------------------------------------*/

class RK3SSP
:
    public timeIntegrator
{

public:

    //- Runtime type information
    TypeName("RK3SSP");

    // Constructor
    RK3SSP(const fvMesh& mesh);


    //- Destructor
    virtual ~RK3SSP();


    // Member Functions

        //- Set ode fields
        virtual void setODEFields(integrationSystem& system);

        //- Update
        virtual void integrate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace timeIntegrators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
