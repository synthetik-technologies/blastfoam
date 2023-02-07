/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::accelerationSchemeList

Description
    Base class for acceleration

SourceFiles
    accelerationSchemeList.C

\*---------------------------------------------------------------------------*/

#ifndef accelerationSchemeList_H
#define accelerationSchemeList_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "accelerationSchemeNew.H"
#include "PtrList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class accelerationSchemeList Declaration
\*---------------------------------------------------------------------------*/

class accelerationSchemeList
:
    public PtrList<accelerationScheme>
{
protected:
// Protected data

    //- Reference to the mesh
    const fvMesh& mesh_;

public:

    TypeName("accelerationSchemeList");

    // Constructor
    accelerationSchemeList
    (
        const word& field,
        const fvMesh& mesh,
        const labelList& patches,
        const dictionary& dict
    );

    //- Destructor
    virtual ~accelerationSchemeList();


    // Member functions

        //- Update residuals
        void updateError();

        //- Relax the field
        void relax(const label iter);

        //- Clear residuals
        void clear();

        // Return the absolute error
        scalar error() const;

        //- Return the relative error
        scalar relError() const;

        //- Read parameters
        void read(const dictionary& dict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
