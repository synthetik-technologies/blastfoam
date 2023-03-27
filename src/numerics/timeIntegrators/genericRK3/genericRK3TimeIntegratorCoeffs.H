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
    Foam::timeIntegrators::genericRK3

Description
    Third order Runge-Kutta time integration with user defined coefficient

SourceFiles
    genericRK3TimeIntegratorCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef genericRK3TimeIntegratorCoeffs_H
#define genericRK3TimeIntegratorCoeffs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrator.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class genericRK3 Declaration
\*---------------------------------------------------------------------------*/

class genericRK3
:
    public timeIntegratorCoeffs
{
    // Private data

        //- Coefficient
        const scalar alpha_;

public:

    //- Runtime type information
    TypeName("genericRK3");

    // Constructor
    genericRK3(Istream& is);
    genericRK3(const scalar alpha);


    //- Destructor
    virtual ~genericRK3();

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
