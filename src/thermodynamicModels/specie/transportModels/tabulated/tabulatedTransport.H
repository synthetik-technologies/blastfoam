/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::tabulatedTransport

Description
    Transport properties package using uniformly-spaced tabulated data for
    viscosity and thermal conductivity vs density and temperature.

Usage
    \table
        Property        | Description
        mu              | Dynamic viscosity vs density and temperature table
        kappa           | Thermal conductivity vs density and temperature table
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        muCoeffs
        {
            rho ( ... );
            T ( ... );
            mu ( ... );
        }

        kappa
        {
            rho ( ... );
            T ( ... );
            kappa ( ... );
        }
    }
    \endverbatim

SourceFiles
    tabulatedTransportI.H
    tabulatedTransport.C

See also
    Foam::lookupTable2D

\*---------------------------------------------------------------------------*/

#ifndef tabulatedTransport_H
#define tabulatedTransport_H

#include "lookupTables2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class tabulatedTransport;

template<class Thermo>
Ostream& operator<<(Ostream&, const tabulatedTransport<Thermo>&);


/*---------------------------------------------------------------------------*\
                     Class tabulatedTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedTransport
:
    public Thermo
{

    // Private Data

        //- Dynamic viscosity table [Pa.s]
        scalarLookupTable2D mu_;

        //- Thermal conductivity table [W/m/K]
        scalarLookupTable2D kappa_;


public:

    // Constructors

        //- Construct from dictionary
        tabulatedTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tabulatedTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulated<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tabulatedTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedTransportI.H"

#ifdef NoRepository
    #include "tabulatedTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
