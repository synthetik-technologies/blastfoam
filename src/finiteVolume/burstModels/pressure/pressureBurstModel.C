/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

\*---------------------------------------------------------------------------*/

#include "pressureBurstModel.H"
#include "burstFvPatchFieldBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace burstModels
{
    defineTypeNameAndDebug(pressure, 0);
    addToRunTimeSelectionTable(burstModel, pressure, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModels::pressure::pressure(const dictionary& dict)
:
    burstModel(dict),
    pName_(dict.lookupOrDefault<word>("pName", "p")),
    pRef_(dict.lookupOrDefault<scalar>("pRef", 0.0)),
    pBurst_(dict.lookup<scalar>("pBurst"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModels::pressure::~pressure()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

bool Foam::burstModels::pressure::update
(
    const fvPatch& patch,
    scalarField& intact
) const
{
    if (burst_)
    {
        return false;
    }

    const fvMesh& mesh = patch.boundaryMesh().mesh();
    bool burst = false;

    if (mesh.foundObject<volScalarField>(pName_))
    {
        const volScalarField& p =
            mesh.lookupObject<volScalarField>(pName_);
        const fvPatchField<scalar>& pp =
            p.boundaryField()[patch.index()];
        scalarField deltaP(pp);

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
            deltaP -= pRef_;
        }

        if (partialBurst_)
        {
            forAll(intact, facei)
            {
                if (deltaP[facei] > pBurst_)
                {
                    intact[facei] = 0.0;
                    burst = true;
                }
            }
            burst_ = gMax(intact) < small;
        }
        else
        {
            // Patch has already burst
            burst = gMax(deltaP) > pBurst_;
            intact = !burst;
            burst_ = burst;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Could not find " << pName_ << " field" << nl
            << "You should be using the impulse functionObject with this " << nl
            << "model" << endl
            << abort(FatalError);
    }
    return returnReduce(burst, orOp<bool>());
}


void Foam::burstModels::pressure::writeData(Ostream& os) const
{
    burstModel::writeData(os);
    writeEntry(os, "pName", pName_);
    writeEntry(os, "pBurst", pBurst_);
    writeEntry(os, "pRef", pRef_);
}


// ************************************************************************* //
