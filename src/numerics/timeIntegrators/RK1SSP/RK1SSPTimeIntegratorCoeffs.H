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
    Foam::timeIntegrators::RK1SSP

Description
    First order, RK1SSP integration

SourceFiles
    RK1SSPTimeIntegratorCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef RK1SSPTimeIntegratorCoeffs_H
#define RK1SSPTimeIntegratorCoeffs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegratorCoeffs.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class RK1SSP Declaration
\*---------------------------------------------------------------------------*/

class RK1SSP
:
    public timeIntegratorCoeffs
{

public:

    //- Runtime type information
    TypeName("RK1SSP");

    // Constructor
    RK1SSP(Istream& is);


    //- Destructor
    virtual ~RK1SSP();

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
            return scalar(nSteps_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace timeIntegrators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
