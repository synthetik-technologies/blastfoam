/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2023
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Third order, strong stability preserving Runge-Kutta method

    References:
    \verbatim
        Spiteri, R.J., Ruuth, S.J. (2003).
        A New Class of Optimal Higher-order Strong-Stability-Preserving
        Time Discretization Methods
        SIAM Journal of Numerical Analysis, 40(2), 469-491.
    \endverbatim

SourceFiles
    RK3SSPTimeIntegratorCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef RK3SSPTimeIntegratorCoeffs_H
#define RK3SSPTimeIntegratorCoeffs_H

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
    public timeIntegratorCoeffs
{

public:

    //- Runtime type information
    TypeName("RK3SSP");

    // Constructor
    RK3SSP(Istream& is);


    //- Destructor
    virtual ~RK3SSP();

    // Public member functions

        //- Set coefficients
        virtual void set
        (
            List<List<scalar>>& as,
            List<List<scalar>>& bs,
            const label index
        ) const;

        //- Return max Courant number
        virtual scalar maxCo() const
        {
            return scalar(nSteps_) - 2.0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace timeIntegrators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
