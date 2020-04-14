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
    Foam::timeIntegrators::RK4SSP

Description
    Fourth order, strong stability preserving Runge-Kutta method

    References:
    \verbatim
        Spiteri, R.J., Ruuth, S.J. (2003).
        A New Class of Optimal Higher-order Strong-Stability-Preserving
        Time Discretization Methods
        SIAM Journal of Numerical Analysis, 40(2), 469-491.
    \endverbatim

SourceFiles
    RK4SSPTimeIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef RK4SSPTimeIntegrator_H
#define RK4SSPTimeIntegrator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrator.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class RK4SSP Declaration
\*---------------------------------------------------------------------------*/

class RK4SSP
:
    public timeIntegrator
{

    // Coefficients
    scalar a20 = 649.0/1600.0;
    scalar b20 = -10890423.0/25193600.0;
    scalar a21 = 951.0/1600.0;
    scalar b21 = 5000.0/7873.0;

    scalar a30 = 53989.0/2500000.0;
    scalar b30 = -102261.0/5000000.0;
    scalar a31 = 4806213.0/20000000.0;
    scalar b31 = -5121.0/20000.0;
    scalar a32 = 23619.0/32000.0;
    scalar b32 = 7873.0/10000.0;

    scalar a40 = 0.2;
    scalar b40 = 0.1;
    scalar a41 = 6127.0/30000.0;
    scalar b41 = 1.0/6.0;
    scalar a42 = 7873.0/30000.0;
    scalar b42 = 0.0;
    scalar a43 = 1.0/3.0;
    scalar b43 = 1.0/6.0;

public:

    //- Runtime type information
    TypeName("RK4SSP");

    // Constructor
    RK4SSP(const fvMesh& mesh);


    //- Destructor
    virtual ~RK4SSP();


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
