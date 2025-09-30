/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2023
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
    Foam::timeIntegrators::RKCK45

Description
    Fifth order Cash-Karp method

    References:
    \verbatim
        Erwin Fehlberg (1969).
        Low-order classical Runge-Kutta formulas with step size control and
        their application to some heat transfer problems.
        NASA Technical Report 315.
    \endverbatim

SourceFiles
    RKCK45TimeIntegratorCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef RKCK45TimeIntegratorCoeffs_H
#define RKCK45TimeIntegratorCoeffs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegratorCoeffs.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class RKCK45 Declaration
\*---------------------------------------------------------------------------*/

class RKCK45
:
    public timeIntegratorCoeffs
{
public:

    //- Runtime type information
    TypeName("RKCK45");

    // Constructor
    RKCK45(Istream& is);


    //- Destructor
    virtual ~RKCK45();

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
