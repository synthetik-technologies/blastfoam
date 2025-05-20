/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::fvMeshTopoChangers::blastList

Description
    Sequence of mesh topology changes applied in order.

SourceFiles
    blastListFvMeshTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef blastListFvMeshTopoChanger_H
#define blastListFvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "PtrDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                   Class fvMeshTopoChangers::blastList Declaration
\*---------------------------------------------------------------------------*/

class blastList
:
    public fvMeshTopoChanger
{
    // Private Data

        PtrDictionary<fvMeshTopoChanger> list_;


public:

    //- Runtime type information
    TypeName("blastList");


    // Constructors

        //- Construct from fvMesh and dictionary
        blastList(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        blastList(const blastList&) = delete;


    //- Destructor
    virtual ~blastList();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Write the mover state
        virtual bool write(const bool write = true) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const list&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
