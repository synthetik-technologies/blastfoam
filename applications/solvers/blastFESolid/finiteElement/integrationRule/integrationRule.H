/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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
    Foam::integrationRule

Description
    List of integration points

SourceFiles
    integrationRuleI.H
    integrationRule.C

\*---------------------------------------------------------------------------*/
#ifndef integrationRule_H
#define integrationRule_H


#include "integrationPoint.H"
#include "ElementType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class integrationRule
:
    public List<integrationPoint>
{
public:

    // List constructors
    using List<integrationPoint>::List;


    //- Stored integration rules
    static FixedList<List<integrationRule>, 8> integrationRules;


    // Member functions

        //- Set mid point of a triangle
        inline void setTriMid
        (
            const label start,
            const scalar w
        );

        //- Set triangle points
        inline void setTri3
        (
            const label start,
            const scalar a,
            const scalar b,
            const scalar w
        );

        inline void setTri3
        (
            const label start,
            const scalar a,
            const scalar w
        );

        inline void setTri3b
        (
            const label start,
            const scalar b,
            const scalar w
        );

        inline void setTri3R
        (
            const label start,
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar w
        );

        inline void setTri3R
        (
            const label start,
            const scalar a,
            const scalar b,
            const scalar w
        );

        inline void setTri6
        (
            const label start,
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar w
        );

        inline void setTri6
        (
            const label start,
            const scalar a,
            const scalar b,
            const scalar w
        );


        //- Set tet points
        inline void setTet3
        (
            const label start,
            const scalar a,
            const scalar b,
            const scalar w
        );

        inline void setTet6
        (
            const label start,
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar w
        );

        inline void setTet6
        (
            const label start,
            const scalar a,
            const scalar w
        );

        inline void setTetMid(const label start, const scalar w);

        inline void setTet4
        (
            const label start,
            const scalar a,
            const scalar w
        );

        inline void setTet4b
        (
            const label start,
            const scalar b,
            const scalar w
        );

        inline void setTet12
        (
            const label start,
            const scalar a,
            const scalar bc,
            const scalar w
        );

        inline void setTet12bc
        (
            const label start,
            const scalar b,
            const scalar c,
            const scalar w
        );

    // Static functions

        //- Return rule for an element type and order
        static const integrationRule& getRule
        (
            const ElementType::Type et,
            const label order
        );

        //- Return point rule for a given order
        static const integrationRule& ptRule(const label order);

        //- Return segment rule for a given order
        static const integrationRule& segRule(const label order);

        //- Return triangle rule for a given order
        static const integrationRule& triRule(const label order);

        //- Return quadrilateral rule for a given order
        static const integrationRule& quadRule(const label order);

        //- Return tetrahedron rule for a given order
        static const integrationRule& tetRule(const label order);

        //- Return hexahedron rule for a given order
        static const integrationRule& hexRule(const label order);

        //- Return prism rule for a given order
        static const integrationRule& prismRule(const label order);

        //- Return pyramid rule for a given order
        static const integrationRule& pyrRule(const label order);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "integrationRuleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
