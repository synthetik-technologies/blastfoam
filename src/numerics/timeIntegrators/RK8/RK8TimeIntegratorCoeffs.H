/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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
    Foam::timeIntegrators::RK8

Description
    Eigth order Runge-Kutta time integration

SourceFiles
    RK8TimeIntegratorCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef RK8TimeIntegratorCoeffs_H
#define RK8TimeIntegratorCoeffs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegratorCoeffs.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class RK8 Declaration
\*---------------------------------------------------------------------------*/

class RK8
:
    public timeIntegratorCoeffs
{

public:

    //- Runtime type information
    TypeName("RK8");

    // Constructor
    RK8(Istream& is);


    //- Destructor
    virtual ~RK8();

    // Public member functions

        //- Set coefficients
        virtual void set
        (
            List<List<scalar>>& as,
            List<List<scalar>>& bs,
            const label index
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace timeIntegrators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
