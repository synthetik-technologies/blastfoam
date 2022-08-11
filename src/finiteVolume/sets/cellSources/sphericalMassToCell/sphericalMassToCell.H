/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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
    Foam::sphericalMassToCell

Description
    Set a sphereCellSet given a mass

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

    // Protected Member Functions

        //- Return the radius
        scalar R(const dictionary& dict) const
        {
            if (mustRead() || volume() > small)
            {
                return
                    cbrt
                    (
                        volume()*3.0/4.0/Foam::constant::mathematical::pi
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
