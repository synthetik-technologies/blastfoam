/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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
    Foam::annulusCylindericalMassToCell

Description
    Set an annulus cylinder given a mass, mean radius, and points

SourceFiles
    annulusCylindericalMassToCell.C

\*---------------------------------------------------------------------------*/

#ifndef annulusCylindericalMassToCell_H
#define annulusCylindericalMassToCell_H

#include "cylinderAnnulusToCell.H"
#include "massToCell.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class annulusCylindericalMassToCell Declaration
\*---------------------------------------------------------------------------*/

class annulusCylindericalMassToCell
:
    public massToCell,
    public cylinderAnnulusToCell
{
// Private data

    //- Mean radius
    scalar radius_;

    //- Span of the solid (rInner = r - deltaR, rOuter = r + deltaR)
    scalar delta_;

    //- Centre point
    point point1_;
    point point2_;

    // Private Member Functions

        void read(const dictionary& dict)
        {
            if (read_)
            {
                return;
            }

            if (mustRead() || volume() > small)
            {
                radius_ = dict.lookup<scalar>("radius");

                if (dict.found("L"))
                {
                    point centre = dict.lookup<point>("centre");
                    vector L = dict.lookup<vector>("L");
                    point1_ = centre - 0.5*L;
                    point2_ = centre + 0.5*L;
                }
                else if (dict.found("point1"))
                {
                    point1_ = dict.lookup<point>("point1");
                    point2_ = dict.lookup<point>("point2");
                }

                scalar L = mag(point1_ - point2_);

                delta_ =
                    volume()
                   /(4.0*L*radius_*Foam::constant::mathematical::pi);
            }
            else
            {
                scalar rInner = dict.lookup<scalar>("innerRadius");
                scalar rOuter = dict.lookup<scalar>("outerRadius");
                point1_ = dict.lookup<point>("point1");
                point2_ = dict.lookup<point>("point2");

                radius_ = 0.5*(rInner + rOuter);
                delta_ = 0.5*(rOuter - rInner);
            }
        }

        //- Return the radius
        vector p1(const dictionary& dict)
        {
            read(dict);

            return point1_;
        }

        vector p2(const dictionary& dict)
        {
            read(dict);

            return point2_;
        }

        scalar outerR(const dictionary& dict)
        {
            read(dict);
            return radius_ + delta_;
        }

        scalar innerR(const dictionary& dict)
        {
            read(dict);
            return radius_ - delta_;
        }

public:

    //- Runtime type information
    TypeName("annulusCylindericalMassToCell");


    // Constructors

        //- Construct from dictionary
        annulusCylindericalMassToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~annulusCylindericalMassToCell();

    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
