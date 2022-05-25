/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::tabulatedIsoSolidTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    tabulatedIsoSolidTransportI.H
    tabulatedIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedIsoSolidTransport_H
#define tabulatedIsoSolidTransport_H

#include "lookupTables1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class tabulatedIsoSolidTransport;

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const tabulatedIsoSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                   Class tabulatedIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedIsoSolidTransport
:
    public Thermo
{
    // Private Data

        //- Lookup table for thermal conductivity
        scalarLookupTable1D kappaTable_;


public:

    // Constructors

        //- Construct from dictionary
        tabulatedIsoSolidTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tabulatedIsoSolidTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedIsoSolidTransport> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulatedIso<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Isotropic thermal conductivity [W/m/K]
        inline scalar kappa
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Un-isotropic thermal conductivity [W/m/K]
        inline vector Kappa
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu
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
            const tabulatedIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedIsoSolidTransportI.H"

#ifdef NoRepository
    #include "tabulatedIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
