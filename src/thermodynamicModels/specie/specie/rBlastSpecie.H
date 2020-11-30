/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Renamed class
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
    Foam::rBlastSpecie

Description

SourceFiles
    rBlastSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef rBlastSpecie_H
#define rBlastSpecie_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "word.H"
#include "scalar.H"
#include "dictionary.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class rBlastSpecie;

inline rBlastSpecie operator+(const rBlastSpecie&, const rBlastSpecie&);
inline rBlastSpecie operator*(const scalar, const rBlastSpecie&);
inline rBlastSpecie operator==(const rBlastSpecie&, const rBlastSpecie&);

/*---------------------------------------------------------------------------*\
                           Class rBlastSpecie Declaration
\*---------------------------------------------------------------------------*/

class rBlastSpecie
{
protected:

    // Protected data

        //- Name of specie
        word specieName_;

        //- Number of moles of this component in the mixture
        scalar Y_;

        //- Molecular weight
        scalar molWeight_;


public:

    //- Runtime type information
    ClassName("specie");



    // Constructors

        //- Construct from components without name
        inline rBlastSpecie(const scalar Y, const scalar molWeight);

        //- Construct from components with name
        inline rBlastSpecie
        (
            const word& name,
            const scalar Y,
            const scalar molWeight
        );

        //- Construct as named copy
        inline rBlastSpecie(const word& name, const rBlastSpecie&);

        //- Construct from dictionary
        rBlastSpecie(const dictionary& dict);


    // Member Functions

        // Access

            //- Name
            inline const word& specieName() const;

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- No of moles of this species in mixture
            inline scalar Y() const;

            //- Gas constant [J/(kg K)]
            inline scalar R() const;


    // Member operators

        inline void operator=(const rBlastSpecie&);
        inline void operator+=(const rBlastSpecie&);
        inline void operator*=(const scalar);


    // Friend operators

        inline friend rBlastSpecie operator+
        (
            const rBlastSpecie&,
            const rBlastSpecie&
        );
        inline friend rBlastSpecie operator*(const scalar, const rBlastSpecie&);
        inline friend rBlastSpecie operator==
        (
            const rBlastSpecie&,
            const rBlastSpecie&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rBlastSpecieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
