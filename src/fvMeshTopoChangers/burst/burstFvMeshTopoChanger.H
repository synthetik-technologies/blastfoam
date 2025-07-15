/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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
    Foam::fvMeshTopoChangers::burst

Description
    TopoChanger that handles the switching of patches when a burst condition is
    met. Either coupled or uncoupled boundaries can be used and how the burst
     criteria is determined is handled by the burstModel.

     Mapping of "new" faces is handles for all

Example
    topoChanger
    {
        type            burst;
        patches
        (
            bag
            {
                intactMasterPatch   bag_master;
                intactSlavePatch    bag_slave;

                burstMasterPatch    bagBurst_master;
                burstSlavePatch     bagBurst_slave;

                burstModel          pressure;
                pBurst              1e5;
                partialBurst        yes;
                useAverage          no;

                patchToPatch        nearest;
            }
        );
    }

SourceFiles
    burstFvMeshTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef burstFvMeshTopoChanger_H
#define burstFvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "polyTopoChange.H"
#include "patchToPatch.H"
#include "burstModel.H"
#include "HashPtrTable.H"
#include "fieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                   Class fvMeshTopoChangers::burst Declaration
\*---------------------------------------------------------------------------*/

class burst
:
    public fvMeshTopoChanger
{
    // Private Data
    class burstData
    {
        //- Is this a coupled set?
        bool coupled_;

        //- Name of intact patches
        const word intactMasterPatchName_;
        const word intactSlavePatchName_;

        //- Names of burst patches
        const word burstMasterPatchName_;
        const word burstSlavePatchName_;

        //- Burst model
        autoPtr<burstModel> burst_;

    public:

        //- Coupled constructor
        burstData
        (
            const word& intactMasterPatchName,
            const word& intactSlavePatchName,

            const word& burstMasterPatchName,
            const word& burstSlavePatchName,

            const dictionary& dict
        )
        :
            coupled_(true),

            intactMasterPatchName_(intactMasterPatchName),
            intactSlavePatchName_(intactSlavePatchName),

            burstMasterPatchName_(burstMasterPatchName),
            burstSlavePatchName_(burstSlavePatchName),

            burst_(burstModel::New(dict, coupled_))
        {}

        //- One-sided constructor
        burstData
        (
            const word& intactPatchName,
            const word& burstPatchName,

            const dictionary& dict
        )
        :
            coupled_(false),

            intactMasterPatchName_(intactPatchName),
            intactSlavePatchName_(word::null),

            burstMasterPatchName_(burstPatchName),
            burstSlavePatchName_(word::null),

            burst_(burstModel::New(dict, coupled_))
        {}

        //- Return is the set is coupled (i.e. master and slave)
        bool coupled() const
        {
            return coupled_;
        }

        //- Return intact master-side patch name
        const word& intactMasterPatch() const
        {
            return intactMasterPatchName_;
        }

        //- Return intact slave-side patch name
        const word& intactSlavePatch() const
        {
            return intactSlavePatchName_;
        }

        //- Return intact patch name
        const word& intactPatch() const
        {
            return intactMasterPatchName_;
        }

        //- Return burst master-side patch name
        const word& burstMasterPatch() const
        {
            return burstMasterPatchName_;
        }

        //- Return burst slave-side patch name
        const word& burstSlavePatch() const
        {
            return burstSlavePatchName_;
        }

        //- Return burst patch name
        const word& burstPatch() const
        {
            return burstMasterPatchName_;
        }

        //- Return burst model
        const burstModel& burst() const
        {
            return burst_();
        }

        //- Access burst model
        burstModel& burst()
        {
            return burst_();
        }
    };

    //- List of burst sets
    PtrList<burstData> patchData_;

    // Private member functions

        //- Store mapped vol boundary fields before mesh changes
        template<class Type>
        void storeVolBoundaries
        (
            const labelList& boundaryMap,
            HashPtrTable<typename VolField<Type>::Boundary>& bfields
        ) const;

        //- Store mapped surface boundary fields before mesh changes
        template<class Type>
        void storeSurfaceBoundaries
        (
            const labelList& boundaryMap,
            HashPtrTable<typename SurfaceField<Type>::Boundary>& bfields
        ) const;

        //- Map vol boundary fields
        template<class Type>
        void mapVolBoundaries
        (
            const List<labelList>& addrs,
            const HashPtrTable<typename VolField<Type>::Boundary>& bfields
        ) const;

        //- Map surface boundary fields
        template<class Type>
        void mapSurfaceBoundaries
        (
            const List<labelList>& addrs,
            const HashPtrTable<typename SurfaceField<Type>::Boundary>& bfields
        ) const;


public:

    //- Runtime type information
    TypeName("burst");


    // Constructors

        //- Construct from fvMesh and dictionary
        burst(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        burst(const burst&) = delete;


    //- Destructor
    virtual ~burst();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const burst&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "burstFvMeshTopoChangerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
