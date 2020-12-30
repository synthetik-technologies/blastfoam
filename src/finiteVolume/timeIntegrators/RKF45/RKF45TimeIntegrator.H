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

    scalar b10 = 0.25;

    scalar b20 = 3.0/32.0;
    scalar b21 = 9.0/32.0;

    scalar b30 = 1932.0/2197.0;
    scalar b31 = -7200.0/2197.0;
    scalar b32 = 7296.0/2197.0;

    scalar b40 = 439.0/216.0;
    scalar b41 = -8.0;
    scalar b42 = 3680.0/513.0;
    scalar b43 = -845.0/4104.0;

    scalar b50 = -8.0/27.0;
    scalar b51 = 2.0;
    scalar b52 = -3544.0/2565.0;
    scalar b53 = 1859.0/4104.0;
    scalar b54 = -11.0/40.0;

    //- Fourth order coefficients
//     scalar b60 = 25.0/216.0;
//     scalar b61 = 0.0;
//     scalar b62 = 1408.0/2565.0;
//     scalar b63 = 2197.0/4104.0;
//     scalar b64 = -1.0/5.0;
//     scalar b65 = 0.0;

    //- Fifth order coefficients
    scalar b60 = 16.0/135.0;
    scalar b61 = 0.0;
    scalar b62 = 6656.0/12825.0;
    scalar b63 = 28561.0/56430.0;
    scalar b64 = -9.0/50.0;
    scalar b65 = 2.0/55.0;


public:

    //- Runtime type information
    TypeName("RKF45");

    // Constructor
    RKF45(const fvMesh& mesh, const label nSteps);


    //- Destructor
    virtual ~RKF45();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace timeIntegrators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
