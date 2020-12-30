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
    Foam::timeIntegrators::RK2SSP

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
    RK2SSP.C

\*---------------------------------------------------------------------------*/

#ifndef RK2SSP_H
#define RK2SSP_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrator.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class RK2SSP Declaration
\*---------------------------------------------------------------------------*/

class RK2SSP
:
    public timeIntegrator
{

public:

    //- Runtime type information
    TypeName("RK2SSP");

    // Constructor
    RK2SSP(phaseCompressibleSystem& fluid);


    //- Destructor
    virtual ~RK2SSP();


    // Member Functions

        //- Update
        virtual void integrate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace timeIntegrators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
