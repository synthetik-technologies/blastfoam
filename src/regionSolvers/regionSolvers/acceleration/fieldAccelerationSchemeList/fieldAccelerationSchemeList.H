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
    Foam::fieldAccelerationSchemeList

Description
    Base class for acceleration

SourceFiles
    fieldAccelerationSchemeList.C

\*---------------------------------------------------------------------------*/

#ifndef fieldAccelerationSchemeList_H
#define fieldAccelerationSchemeList_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldAccelerationScheme.H"
#include "regionSolverList.H"
#include "PtrListDictionary.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fieldAccelerationSchemeList Declaration
\*---------------------------------------------------------------------------*/

class fieldAccelerationSchemeList
:
    public PtrListDictionary<fieldAccelerationScheme>
{
protected:
// Protected data

    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Reference to the regionSolverList object
    const regionSolverList& regions_;

public:

    TypeName("fieldAccelerationSchemeList");

    // Constructor
    fieldAccelerationSchemeList
    (
        const fvMesh& mesh,
        const regionSolverList& regions
    );

    //- Destructor
    virtual ~fieldAccelerationSchemeList();


    // Member functions

        //- Add a field for convergence
        void addField
        (
            const word& fieldName,
            const labelList& patches
        );

        //- Add a field for convergence
        void addField(const word& fieldName);

        //- Update residuals
        void updateError();

        //- Relax the field
        void relax(const label iter);

        //- Clear residuals
        void clear();

        //- Return convergence of all fields
        Convergence convergence() const;

        //- Read parameters
        void read(const dictionary& dict);

        //- Print current errors
        void print(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
