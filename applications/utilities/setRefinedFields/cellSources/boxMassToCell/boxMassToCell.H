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
    Foam::boxMassToCell

Description
    A topoSetSource to select cells based on cell centres inside a list
    of cylinders.

SourceFiles
    boxMassToCell.C

\*---------------------------------------------------------------------------*/

#ifndef boxMassToCell_H
#define boxMassToCell_H

#include "boxToCell.H"
#include "massToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class boxMassToCell Declaration
\*---------------------------------------------------------------------------*/

class boxMassToCell
:
    public massToCell,
    public boxToCell
{

    //- Return the bound box
    treeBoundBoxList bounds(const dictionary& dict) const
    {
        vector p1, p2;
        if (mustRead() || volume() > small)
        {
            vector aspectRatio(dict.lookup<vector>("aspectRatio"));
            scalar x =
                cbrt
                (
                    volume()
                   /(aspectRatio.x()*aspectRatio.y()*aspectRatio.z())
                );

            vector span(x*aspectRatio);
            p1 = centre() - span/2.0;
            p2 = centre() + span/2.0;
            scalar v = span.x()*span.y()*span.z();
        }
        else
        {
            vector span(dict.lookup<vector>("span"));
            p1 = centre() - span/2.0;
            p2 = centre() - span/2.0;
        }

        return treeBoundBoxList({treeBoundBox(p1, p2)});
    }


public:

    //- Runtime type information
    TypeName("boxMassToCell");


    // Constructors

        //- Construct from dictionary
        boxMassToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~boxMassToCell();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
