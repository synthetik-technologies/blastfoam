/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "burstFvPatchFieldBase.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstFvPatchFieldBase, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstFvPatchFieldBase::burstFvPatchFieldBase
(
    const fvPatch& p
)
:
    burstFvPatch_
    (
        const_cast<burstFvPatchBase&>
        (
            dynamicCast<const burstFvPatchBase>(p)
        )
    ),
    burstPolyPatch_(burstFvPatch_.burstPolyPatch()),

    // Default to unblocked (i.e. for temporary fields)
    unblock_(false),
    block_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstFvPatchFieldBase::~burstFvPatchFieldBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstFvPatchFieldBase::update()
{
    if (!burstPolyPatch_.uptoDate(true))
    {
        const polyPatch& p = burstPolyPatch_.patch();
        const polyMesh& mesh = p.boundaryMesh().mesh();

        scalarField deltaP(p.size(), 0.0);
        scalarField deltaImp(p.size(), 0.0);
        bool needUpdate = false;

        if (burstPolyPatch_.usePressure())
        {
            if
            (
                mesh.template
                foundObject<volScalarField>(burstPolyPatch_.pName())
            )
            {
                const volScalarField& pF =
                    mesh.template
                    lookupObject<volScalarField>(burstPolyPatch_.pName());
                const fvPatchField<scalar>& pp =
                    pF.boundaryField()[p.index()];

                if (pp.coupled())
                {
                    // Unblock the patch so that the coupled
                    // boundary type is used to find the neighbour field
                    if (isA<burstFvPatchFieldBase>(pp))
                    {
                        dynamicCast<const burstFvPatchFieldBase>
                        (
                            pp
                        ).unblock(true);
                    }
                    deltaP =
                        mag
                        (
                            pp.patchInternalField()
                          - pp.patchNeighbourField()
                        );
                    if (isA<burstFvPatchFieldBase>(pp))
                    {
                        // Reset the unblock flag
                        dynamicCast<const burstFvPatchFieldBase>
                        (
                            pp
                        ).unblock(false);
                    }
                }
                else
                {
                    deltaP = mag(pp - burstPolyPatch_.pRef());
                }
                needUpdate = true;
            }
//             else
//             {
//                 WarningInFunction
//                     << "Could not find " << burstPolyPatch_.pName()
//                     << ", neglecting pressure. " << endl;
//             }
        }
        if (burstPolyPatch_.useImpulse())
        {
            if
            (
                mesh.template
                foundObject<volScalarField>(burstPolyPatch_.impulseName())
            )
            {
                const volScalarField& impF =
                    mesh.template
                    lookupObject<volScalarField>(burstPolyPatch_.impulseName());
                const fvPatchField<scalar>& pimp =
                    impF.boundaryField()[p.index()];
                if (pimp.coupled())
                {
                    // Unblock the patch so that the coupled
                    // boundary type is used to find the neighbour field
                    if (isA<burstFvPatchFieldBase>(pimp))
                    {
                        dynamicCast<const burstFvPatchFieldBase>
                        (
                            pimp
                        ).unblock(true);
                    }
                    deltaImp =
                        mag
                        (
                            pimp.patchInternalField()
                          - pimp.patchNeighbourField()
                        );

                    // Reset the unblock flag
                    if (isA<burstFvPatchFieldBase>(pimp))
                    {
                        dynamicCast<const burstFvPatchFieldBase>
                        (
                            pimp
                        ).unblock(false);
                    }
                }
                else
                {
                    deltaImp = mag(pimp);
                }
                needUpdate = true;
            }
//             else
//             {
//                 WarningInFunction
//                     << "Could not find " << burstPolyPatch_.impulseName()
//                     << ", neglecting impulse. " << endl;
//             }
        }

        if (returnReduce(needUpdate, orOp<bool>()))
        {
            burstFvPatch_.update(deltaP, deltaImp);
        }
    }
}

// ************************************************************************* //
