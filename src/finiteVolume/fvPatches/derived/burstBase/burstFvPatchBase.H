/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::burstFvPatchBase

Description
    Base class for burst fvPatches

SourceFiles
    burstFvPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef burstFvPatchBase_H
#define burstFvPatchBase_H

#include "cyclicFvPatch.H"
#include "burstPolyPatchBase.H"
#include "fvBoundaryMesh.H"
#include "fvPatchFieldMapper.H"
#include "regIOobject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class burstFvPatchParent Declaration
\*---------------------------------------------------------------------------*/
class burstFvPatchParent
:
    public regIOobject
{
// Private data

    //- Constant reference to the mesh
    const fvMesh& mesh_;

    //- Do the field need to be read (i.e. are we starting)
    mutable bool needToRead_;

    //- Intact field not owned by this object
    //  the field is not deleted if not owned my this object
    mutable autoPtr<volScalarField> intactPtr_;

    //- Make the intact field
    void makeIntact() const;

public:

    //- Runtime type information
    TypeName("burstFvPatchParent");

    // Constructor
    burstFvPatchParent(const fvMesh& mesh);

    // Destructor
    virtual ~burstFvPatchParent();


    //- Return the intact field
    const volScalarField& intact() const
    {
        if (!intactPtr_.valid())
        {
            makeIntact();
        }
        return intactPtr_();
    }

    //- Return the intact field for a patch
    const tmp<scalarField> intact(const label patchi) const
    {
        if (!intactPtr_.valid())
        {
            makeIntact();

            if (!intactPtr_.valid())
            {
                return tmp<scalarField>
                (
                    new scalarField(mesh_.boundaryMesh()[patchi].size(), 1.0)
                );
            }
        }
        return intactPtr_->boundaryField()[patchi];
    }

    //- Update the intact field
    void update(const label patchi, const scalarField& intact);

    //- Return a reference to the object stored in mesh
    static burstFvPatchParent& New(const fvMesh& mesh);

    virtual bool writeData(Ostream&) const
    {
        return true;
    }
};


/*---------------------------------------------------------------------------*\
                        Class burstFvPatchBase Declaration
\*---------------------------------------------------------------------------*/

class burstFvPatchBase
{
protected:
    // Protected Data

        //- Reference to the fvPatch
        const fvPatch& patch_;

        //- Constant reference to the burstCyclicPolyPatch
        burstPolyPatchBase& burstPolyPatch_;

        //- Reference to the parent
        burstFvPatchParent& parent_;

        //- Current time index
        mutable label curTimeIndex_;

        //- Update the mesh weights and deltas
        void updateDeltas();

        //- Make weights
        virtual void makeWeights(scalarField&) const = 0;


public:

    friend burstFvPatchParent;

    //- Runtime type information
    TypeName("burstFvPatchBase");

    // Constructors

        //- Construct from polyPatch
        burstFvPatchBase(const fvPatch& patch)
        :
            patch_(patch),
            burstPolyPatch_
            (
                const_cast<burstPolyPatchBase&>
                (
                    refCast<const burstPolyPatchBase>(patch.patch()))
            ),
            parent_
            (
                burstFvPatchParent::New(patch_.boundaryMesh().mesh())
            ),
            curTimeIndex_(patch.boundaryMesh().mesh().time().timeIndex())
        {}

        burstFvPatchBase(const burstFvPatchBase&) = delete;

    // Destructor
    virtual ~burstFvPatchBase()
    {}

    // Member Functions

        //- Return the fvPatch
        const fvPatch& patch() const
        {
            return patch_;
        }

        //- Return the burstPolyPatchBase
        const burstPolyPatchBase& burstPolyPatch() const
        {
            return burstPolyPatch_;
        }

        //- Return the intact field
        tmp<scalarField> intact() const
        {
            return parent_.intact(patch_.index());
        }

        //- Update the intact field and update mesh weight and deltas if
        //  a change happens
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
