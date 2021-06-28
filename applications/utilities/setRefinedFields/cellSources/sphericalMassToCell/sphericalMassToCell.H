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
    Foam::sphericalMassToCell

Description
    A topoSetSource to select cells based on cell centres inside sphere.

SourceFiles
    sphericalMassToCell.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalMassToCell_H
#define sphericalMassToCell_H

#include "sphereToCell.H"
#include "massToCell.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sphericalMassToCell Declaration
\*---------------------------------------------------------------------------*/

class sphericalMassToCell
:
    public massToCell,
    public sphereToCell
{

    // Private Member Functions

        //- Return the radius
        scalar R(const dictionary& dict) const
        {
            if (mustRead() || volume() > small)
            {
                return
                    cbrt
                    (
                        volume()*3.0/4.0*Foam::constant::mathematical::pi
                    );
            }
            return dict.lookup<scalar>("radius");
        }


public:

    //- Runtime type information
    TypeName("sphericalMassToCell");


    // Constructors

        //- Construct from dictionary
        sphericalMassToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sphericalMassToCell();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
