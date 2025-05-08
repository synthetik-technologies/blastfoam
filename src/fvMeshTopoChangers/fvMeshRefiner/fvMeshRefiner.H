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
    Foam::fvMeshRefiner

SourceFiles
    fvMeshRefiner.C

Description
    A general class for the selection of mesh refinement method

\*---------------------------------------------------------------------------*/

#ifndef fvMeshRefiner_H
#define fvMeshRefiner_H

#include "errorEstimator.H"
#include "polyMeshRefiner.H"
#include "fvMeshTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class fvMeshRefiner Declaration
\*---------------------------------------------------------------------------*/

class fvMeshRefiner
:
    public fvMeshTopoChanger
{
protected:

        //- Error estimator
        autoPtr<errorEstimator> error_;

        //- polyMesh refiner
        autoPtr<polyMeshRefiner> refiner_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;


public:

    //- Runtime type information
    TypeName("fvMeshRefiner");


    // Constructors

        //- Construct from mesh and dictionary
        fvMeshRefiner
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct from mesh, dictionary, and refiner type
        fvMeshRefiner
        (
            fvMesh& mesh,
            const dictionary& dict,
            const word& refinerType
        );

        //- Disallow default bitwise copy construction
        fvMeshRefiner(const fvMeshRefiner&) = delete;


    //- Destructor
    virtual ~fvMeshRefiner();


    // Member Functions

        // Access

            //- Is the mesh refining
            inline bool isRefining() const
            {
                return refiner_->isRefining();
            }

            //- Is the mesh unrefining
            inline bool isUnrefining() const
            {
                return refiner_->isUnrefining();
            }

            //- Return the polyMeshRefiner
            inline const polyMeshRefiner& refiner() const
            {
                return refiner_();
            }

            //- Access the polyMeshRefiner
            inline polyMeshRefiner& refiner()
            {
                return refiner_();
            }

            //- Return the number of buffer layers
            inline label nRefinementBufferLayers() const
            {
                return refiner_->nRefinementBufferLayers_;
            }

            //- Return the number of buffer layers
            inline label nUnrefinementBufferLayers() const
            {
                return refiner_->nUnrefinementBufferLayers_;
            }

            //- Force refinement/unrefinement
            inline void setForce(const bool force)
            {
                refiner_->setForce(force);
            }

            //- Return the cell level
            inline const labelList& cellLevel() const
            {
                return refiner_->cellLevel();
            }

            //- Return the point level
            inline const labelList& pointLevel() const
            {
                return refiner_->pointLevel();
            }

            //- Return the location mapper
            inline const locationMapper& locMapper() const
            {
                return locationMapper::New(mesh());
            }

            //- Return the errorEstimator
            const errorEstimator& error() const
            {
                return error_();
            }

            //- Access the errorEstimator
            errorEstimator& error()
            {
                return error_();
            }

            //- Tell the location mapper to store data
            inline void needLocationMap() const
            {
                locMapper().needMap();
            }

        //- Refine the mesh
        virtual bool update();

        //- Refine the mesh
        virtual bool refine
        (
            const scalarField& error,
            const labelList& maxCellLevel = labelList(),
            const scalar lowerRefineLevel = sqrt(small),
            const scalar upperRefineLevel = great,
            const scalar unrefineLevel = -sqrt(small)
        );

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Writing

        //- Write
        virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshRefiner&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
