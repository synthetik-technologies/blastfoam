/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |   Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
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
    Foam::adaptiveFvMesh

SourceFiles
    adaptiveFvMesh.C

Authors
    T.G. Voskuilen ( https://github.com/tgvoskuilen/meshBalancing )
    Daniel Deising <deising@mma.tu-darmstadt.de>
    Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
    All rights reserved.

Description
    A fvMesh with run-time load balancing.

    Updated to OpenFOAM-7 from the github version at:
    https://github.com/ElsevierSoftwareX/SOFTX_2018_143.git

    Reference:
    Rettenmaier, D., Deising, D., Ouedraogo, Y., Gjonaj, E., De Gersem, H., Bothe, D., Tropea, C., Marschall, H.
    "Load balanced 2D and 3D adaptive mesh refinement in OpenFOAM" (2019),
    SoftwareX, volume 10, 2352-7110,

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

    Additional modifications not part of the original work include the use of
    error estimators, improved stability with castellated mesh, and fewer
    required user inputs.

\*---------------------------------------------------------------------------*/

#ifndef adaptiveFvMesh_H
#define adaptiveFvMesh_H

#include "dynamicBlastFvMesh.H"
#include "errorEstimator.H"
#include "fvMeshRefiner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class adaptiveFvMesh Declaration
\*---------------------------------------------------------------------------*/

class adaptiveFvMesh
:
    public dynamicBlastFvMesh
{
protected:

        autoPtr<errorEstimator> error_;

        //- Mesh cutting engine
        autoPtr<fvMeshRefiner> refiner_;

        //- Fluxes to map
        HashTable<word> correctFluxes_;


        // Protected member functions

            //- Read from dynamicFvMesh
            void readDict();

            //- Overload update mesh to include other methods
            virtual void updateMesh(const mapPolyMesh& mpm);

            //- Distribute additional information
            virtual void distribute(const mapDistributePolyMesh& map);

            //- Map non-flux surface<Type>Fields for new internal faces
            //  (from cell splitting)
            template <class T>
            void mapNewInternalFaces(const labelList& faceMap)
            {}


public:

    //- Runtime type information
    TypeName("adaptiveFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit adaptiveFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        adaptiveFvMesh(const adaptiveFvMesh&) = delete;


    //- Destructor
    virtual ~adaptiveFvMesh();


    // Member Functions

        //- Map all fields in time using given map
        virtual void mapFields(const mapPolyMesh& mpm);

        //- Direct access to the refinement engine
        const fvMeshRefiner& refiner() const
        {
            return refiner_();
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Refine the mesh
        virtual bool refine();


    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write = true
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const adaptiveFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
