/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::cylindericalMassToCell

Description
    A topoSetSource to select cells based on cell centres inside a
    cylinder defined by a mass and an length to diameter ratio.

SourceFiles
    cylindericalMassToCell.C

\*---------------------------------------------------------------------------*/

#ifndef cylindericalMassToCell_H
#define cylindericalMassToCell_H

#include "cylinderToCell.H"
#include "massToCell.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cylindericalMassToCell Declaration
\*---------------------------------------------------------------------------*/

class cylindericalMassToCell
:
    public massToCell,
    public cylinderToCell
{
    // Private data
    scalar r_;
    scalar L_;
    vector dir_;

    // Private Member Functions

        void read(const dictionary& dict)
        {
            if (read_)
            {
                return;
            }

            if (mustRead() || volume() > small)
            {
                scalar LbyD(dict.lookup<scalar>("LbyD"));
                vector dir(dict.lookup<vector>("direction"));
                r_ =
                    cbrt
                    (
                        volume()/(2.0*LbyD*Foam::constant::mathematical::pi)
                    );
                L_ = 2.0*r_*LbyD;
                dir_ = dir;
            }
            else
            {
                r_ = dict.lookup<scalar>("radius");
                vector L(dict.lookup<vector>("L"));
                L_ = mag(L);
                dir_ = cmptMag(L)/L_;
            }
        }

        //- Return the radius
        vector p1(const dictionary& dict)
        {
            read(dict);

            return centre() - dir_*L_/2.0;
        }

        vector p2(const dictionary& dict)
        {
            read(dict);

            return centre() + dir_*L_/2.0;
        }

        scalar R(const dictionary& dict)
        {
            read(dict);
            return r_;
        }

public:

    //- Runtime type information
    TypeName("cylindericalMassToCell");


    // Constructors

        //- Construct from dictionary
        cylindericalMassToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cylindericalMassToCell();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
