/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
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
    Foam::dynamicRefineBalancedFvMesh

SourceFiles
    dynamicRefineBalancedFvMesh.C
    dynamicRefineBalancedFvMeshTemplates.C

Authors
    T.G. Voskuilen ( https://github.com/tgvoskuilen/meshBalancing )
    Daniel Deising <deising@mma.tu-darmstadt.de>
    Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
    All rights reserved.

Description
    A fvMesh with run-time load balancing.

    You may refer to this software as:
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier (main developer).

    Method Development and Intellectual Property :
        T.G. Voskuilen (Purdue University)
        Timothée Pourpoint <timothee@purdue.edu> (Purdue University)
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
        Daniel Deising <deising@mma.tu-darmstadt.de>
        Holger Marschall <marschall@csi.tu-darmstadt.de>
        Dieter Bothe <bothe@csi.tu-darmstadt.de>
        Cameron Tropea <ctropea@sla.tu-darmstadt.de>

        School of Aeronautics and Astronautics Purdue University
        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineBalancedFvMesh_H
#define dynamicRefineBalancedFvMesh_H


#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "hexRef8.H"
#include "PackedBoolList.H"
#include "Switch.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicRefineBalancedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefineBalancedFvMesh
:
    public dynamicRefineFvMesh
{

private:


        //- Disallow default bitwise copy construct
        dynamicRefineBalancedFvMesh(const dynamicRefineBalancedFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineBalancedFvMesh&);

        //-
        label topParentID(label p);

        //-
        bool rebalance_;

public:

    //- Runtime type information
    TypeName("dynamicRefineBalancedFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefineBalancedFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicRefineBalancedFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Template to update all volField boundaries
        template<class Type> void correctBoundaries();

        //-
        bool rebalance()
        {
            return rebalance_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicRefineBalancedFvMeshTemplates.C"
#endif

#endif

// ************************************************************************* //
