/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::integrationPoint

Description
    Integration point containing vector location and scalar weight

SourceFiles
    integrationPointI.H
    integrationPoint.C

\*---------------------------------------------------------------------------*/

#ifndef integrationPoint_H
#define integrationPoint_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class integrationPoint;

Istream& operator>>(Istream&, integrationPoint&);
Ostream& operator<<(Ostream&, const integrationPoint&);


/*---------------------------------------------------------------------------*\
                    Class integrationPoint Declaration
\*---------------------------------------------------------------------------*/

class integrationPoint
:
    public vector
{
    //- Weight
    scalar w_;


public:


    // Constructors

        //- Null constructor
        inline integrationPoint();

        //- Construct from a point and a weight
        inline integrationPoint(const vector& pt, const scalar w);

        //- Construct from a 3 vector components and a weight
        inline integrationPoint
        (
            const scalar x,
            const scalar y,
            const scalar z,
            const scalar w
        );

        //- Construct from istream
        integrationPoint(Istream& is);

    //- Destructor
    inline ~integrationPoint();


    // Member functions

        // Access
            //- Return weight
            inline scalar w() const;

            //- Access weight
            inline scalar& w();



        //- Set x value only
        inline void set1(const scalar x);

        //- Set x value and weight
        inline void set1(const scalar x, const scalar w);

        //- Set x and y vales
        inline void set2
        (
            const scalar x,
            const scalar y
        );

        //- Set x and y values and weight
        inline void set2
        (
            const scalar x,
            const scalar y,
            const scalar w
        );

        //- Set x, y, and z values
        inline void set3
        (
            const scalar x,
            const scalar y,
            const scalar z
        );

        //- Set x, y, and z values and weight
        inline void set3
        (
            const scalar x,
            const scalar y,
            const scalar z,
            const scalar w
        );


    // IOstream operators

        friend Istream& operator>>(Istream&, integrationPoint&);

        friend Ostream& operator<<(Ostream&, const integrationPoint&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "integrationPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
