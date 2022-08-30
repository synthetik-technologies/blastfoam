/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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
    Foam::fvMeshPolyRefiner

SourceFiles
    fvMeshPolyRefiner.C

Description
    A fvMesh with run-time load balancing.

\*---------------------------------------------------------------------------*/

#ifndef fvMeshPolyRefiner_H
#define fvMeshPolyRefiner_H

#include "fvMeshRefiner.H"
#include "refinement.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fvMeshPolyRefiner Declaration
\*---------------------------------------------------------------------------*/

class fvMeshPolyRefiner
:
    public fvMeshRefiner
{
    //- Refiner
    autoPtr<refinement> refiner_;

    // Private member functions

        //- Select points to unrefine
        labelList selectUnrefinePoints
        (
            const scalar unrefineLevel,
            const PackedBoolList& markedCell,
            const scalarField& vFld
        ) const;


public:

    //- Runtime type information
    TypeName("polyRefiner");


    // Constructors

        //- Construct from fvMesh
        explicit fvMeshPolyRefiner(fvMesh& mesh);

        //- Construct from fvMesh and dictionary
        fvMeshPolyRefiner
        (
            fvMesh& mesh,
            const dictionary& dict,
            const bool force = false,
            const bool read = true
        );

        //- Disallow default bitwise copy construction
        fvMeshPolyRefiner(const fvMeshPolyRefiner&) = delete;


    //- Destructor
    virtual ~fvMeshPolyRefiner();


    // Member Functions

        // Access

            //- Return the cell level
            virtual const labelList& cellLevel() const
            {
                return refiner_->cellLevel();
            }

            //- Return the point level
            virtual const labelList& pointLevel() const
            {
                return refiner_->pointLevel();
            }

            //- Return the refiner
            const refinement& refiner() const
            {
                return refiner_();
            }

            //- Return the location mapper
            virtual const locationMapper& locMapper() const
            {
                return refiner_->locMapper();
            }


        //- Read the projection parameters from dictionary
        void readDict(const dictionary& parentDict);

        //- Refine the mesh
        virtual bool refine
        (
            const scalarField& error,
            const labelList& maxCellLevel = labelList(),
            const scalar lowerRefineLevel = small,
            const scalar upperRefineLevel = great,
            const scalar unrefineLevel = -small
        );

        //- Overload update mesh to include other methods
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Distribute additional information
        virtual void distribute(const mapDistributePolyMesh& map);


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
        void operator=(const fvMeshPolyRefiner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
