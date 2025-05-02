/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::fvMeshTopoChangers::raw

Description
    Prototype fvMeshTopoChanger without any added functionality.

    Can be used as the starting point for any polyTopoChanger based
    topo-changer.

    Note: run without FOAM_SETNAN. Temporary has unitialised patch
    data when faces are created out of nothing.

SourceFiles
    burstFvMeshTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef burstFvMeshTopoChanger_H
#define burstFvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "polyTopoChange.H"
#include "patchToPatch.H"
#include "burstModel.H"

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
        bool coupled_;

        const word intactMasterPatchName_;
        const word intactSlavePatchName_;

        const word burstMasterPatchName_;
        const word burstSlavePatchName_;

        autoPtr<burstModel> burst_;

    public:

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

        bool coupled() const
        {
            return coupled_;
        }

        const word& intactMasterPatch() const
        {
            return intactMasterPatchName_;
        }
        const word& intactSlavePatch() const
        {
            return intactSlavePatchName_;
        }
        const word& intactPatch() const
        {
            return intactMasterPatchName_;
        }

        const word& burstMasterPatch() const
        {
            return burstMasterPatchName_;
        }
        const word& burstSlavePatch() const
        {
            return burstSlavePatchName_;
        }
        const word& burstPatch() const
        {
            return burstMasterPatchName_;
        }

        const burstModel& burst() const
        {
            return burst_();
        }

        burstModel& burst()
        {
            return burst_();
        }
    };

    PtrList<burstData> patchData_;


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

#endif

// ************************************************************************* //
