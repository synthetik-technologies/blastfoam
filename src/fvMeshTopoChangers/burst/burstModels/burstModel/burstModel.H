/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
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
    Foam::burstModel

Description
    Base class for time integration

SourceFiles
    burstModel.C
    burstModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef burstModel_H
#define burstModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "runTimeSelectionTables.H"
#include "patchToPatch.H"
#include "fvPatch.H"
#include "distributionMap.H"
#include "remote.H"


namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class burstModel Declaration
\*---------------------------------------------------------------------------*/

class burstModel
{
protected:
// Protected variables

    //- Constructing dictionary
    const dictionary dict_;

    //- Is this for a coupled set of patches
    word patchToPatchType_;

    //- Is a partial burst used
    bool partialBurst_;

    //- Is the delta across coupled patches used
    bool useDelta_;

    //- Use average pressure to determine burst (only for non partial burst)
    bool useAverage_;

    //- Split patch into connect regions, non-partial burst applied to regions
    bool regionize_;

    //- Log information
    bool log_;

    //- Is an update needed
    mutable bool needUpdate_;
    mutable bool needRegionUpdate_;

    //- Patch to patch storage
    mutable autoPtr<patchToPatch> mappingPtr_;

    //- List of region faces
    mutable List<labelList> patchRegionFaces_;
    mutable List<labelList> patch1RegionFaces_;
    mutable List<labelList> patch2RegionFaces_;

    // Protected member function

        //- Helper function to find faces that have burst
        bool findBurstFaces
        (
            const scalarField& pf1,
            const scalarField& magSf1,
            const scalarField& pf2,
            const scalarField& magSf2,
            const scalar burstValue,
            labelHashSet& master,
            labelHashSet& slave
        ) const;

        bool findBurstFaces
        (
            const scalarField& pf,
            const scalarField& magSf,
            const scalar burstValue,
            labelHashSet& faces
        ) const;

        //- Split single patch in to regions
        void regionize(const polyPatch& patch) const;

        //- Split single patch in to regions
        void regionize
        (
            const polyPatch& patch,
            List<labelList>& patchRegionFaces
        ) const;

        //- Split two coupled patch in to regions
        void regionize
        (
            const polyPatch& patch1,
            const polyPatch& patch2
        ) const;

        //- Find faces that need to be sent to other processors
        labelListList facesToSend
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch
        ) const;

        //- Distribute patch
        List<remote> distributePatch
        (
            const distributionMap& map,
            const primitivePatch& patch,
            autoPtr<PrimitivePatch<faceList, pointField>>& localPatchPtr
        ) const;

public:

    //- Runtime type information
    TypeName("burstModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            burstModel,
            dictionary,
            (
                const dictionary& dict,
                const bool coupled
            ),
            (dict, coupled)
        );

    //- Constructor
    burstModel
    (
        const dictionary& dict,
        const bool coupled
    );


    //- Destructor
    virtual ~burstModel();


    // Selectors

        //- Return a new burstModel object
        static autoPtr<burstModel> New
        (
            const dictionary& dict,
            const bool coupled
        );

        //- Clone
        virtual autoPtr<burstModel> clone() const = 0;


    // Member Functions

        //- Set if an update is needed
        void needUpdate() const;

        //- Update the patch to patch mapping
        void updateMapping
        (
            const polyPatch& patch1,
            const polyPatch& patch2
        ) const;

        //- Find and set face indices to change, return if any added
        virtual bool facesToChange
        (
            const fvPatch& patch1,
            const fvPatch& patch2,
            labelList& masterFaces,
            labelList& slaveFaces
        ) const = 0;

        //- Find and set face indices to change, return if any added
        virtual bool facesToChange
        (
            const fvPatch& patch,
            labelList& faces
        ) const = 0;

        //- Return if burst occurs
        //  1 owner side
        //  -1 neighbour side
        // 0 no burst
        virtual label update
        (
            const objectRegistry& obr,
            const labelList& ownCells,
            const labelList& neiCells,
            const scalarField& V
        ) const = 0;

        //- Write
        virtual void writeData(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
