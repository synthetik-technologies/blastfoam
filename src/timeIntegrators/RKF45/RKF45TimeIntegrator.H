/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
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
    Foam::timeIntegrators::RKF45
Description
    Fifth order Runge–Kutta–Fehlberg method
    References:
    \verbatim
        Erwin Fehlberg (1969).
        Low-order classical Runge-Kutta formulas with step size control and
        their application to some heat transfer problems.
        NASA Technical Report 315.
    \endverbatim
SourceFiles
    RKF45TimeIntegrator.C
\*---------------------------------------------------------------------------*/

#ifndef RKF45TimeIntegrator_H
#define RKF45TimeIntegrator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrator.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class RKF45 Declaration
\*---------------------------------------------------------------------------*/

class RKF45
:
    public timeIntegrator
{

    // Coefficients

    scalar a10 = 0.25;

    scalar a20 = 3.0/32.0;
    scalar a21 = 9.0/32.0;

    scalar a30 = 1932.0/2197.0;
    scalar a31 = -7200.0/2197.0;
    scalar a32 = 7296.0/2197.0;

    scalar a40 = 439.0/216.0;
    scalar a41 = -8.0;
    scalar a42 = 3680.0/513.0;
    scalar a43 = -845.0/4104.0;

    scalar a50 = -8.0/27.0;
    scalar a51 = 2.0;
    scalar a52 = -3544.0/2565.0;
    scalar a53 = 1859.0/4104.0;
    scalar a54 = -11.0/40.0;

    //- Fourth order coefficients
    scalar a604 = 25.0/216.0;
    scalar a614 = 0.0;
    scalar a624 = 1408.0/2565.0;
    scalar a634 = 2197.0/4104.0;
    scalar a644 = -1.0/5.0;
    scalar a654 = 0.0;

    //- Fifth order coefficients
    scalar a605 = 16.0/135.0;
    scalar a615 = 0.0;
    scalar a625 = 6656.0/12825.0;
    scalar a635 = 28561.0/56430.0;
    scalar a645 = -9.0/50.0;
    scalar a655 = 2.0/55.0;


public:

    //- Runtime type information
    TypeName("RKF45");

    // Constructor
    RKF45(const fvMesh& mesh);


    //- Destructor
    virtual ~RKF45();


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
