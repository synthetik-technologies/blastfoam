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
    Foam::timeIntegrators::RK4

Description
    Fourth order Runge-Kutta method

SourceFiles
    RK4TimeIntegratorCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef RK4TimeIntegratorCoeffs_H
#define RK4TimeIntegratorCoeffs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegratorCoeffs.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class RK4 Declaration
\*---------------------------------------------------------------------------*/

class RK4
:
    public timeIntegratorCoeffs
{

public:

    //- Runtime type information
    TypeName("RK4");

    // Constructor
    RK4(Istream& is);


    //- Destructor
    virtual ~RK4();

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
