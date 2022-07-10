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

#include "pressureImpulseBurstModel.H"
#include "burstFvPatchFieldBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace burstModels
{
    defineTypeNameAndDebug(pressureImpulse, 0);
    addToRunTimeSelectionTable(burstModel, pressureImpulse, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModels::pressureImpulse::pressureImpulse
(
    const dictionary& dict
)
:
    pressure(dict),
    impulseName_(dict.lookupOrDefault<word>("impulseName", "impulse")),
    impulseBurst_(dict.lookup<scalar>("impulseBurst"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModels::pressureImpulse::~pressureImpulse()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

bool Foam::burstModels::pressureImpulse::update
(
    const fvPatch& patch,
    scalarField& intact
) const
{
    bool burst = pressure::update(patch, intact);
    if (burst_)
    {
        return burst;
    }
    const fvMesh& mesh = patch.boundaryMesh().mesh();
    if (mesh.foundObject<volScalarField>(impulseName_))
    {
        const volScalarField& impulse =
            mesh.lookupObject<volScalarField>(impulseName_);
        const fvPatchField<scalar>& pimp =
            impulse.boundaryField()[patch.index()];
        scalarField deltaImp(this->patchField(pimp));

        if (partialBurst_)
        {
            forAll(intact, facei)
            {
                if (intact[facei] > 0.5 && deltaImp[facei] > impulseBurst_)
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
            burst = gMax(deltaImp) > impulseBurst_;
            intact = !burst;
            burst_ = burst;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Could not find " << impulseName_ << " field" << nl
            << "You should be using the impulse functionObject with this " << nl
            << "model" << endl
            << abort(FatalError);
    }
    return returnReduce(burst, orOp<bool>());
}


void Foam::burstModels::pressureImpulse::writeData(Ostream& os) const
{
    pressure::writeData(os);
    writeEntry(os, "impulseName", impulseName_);
    writeEntry(os, "impulseBurst", impulseBurst_);
}


// ************************************************************************* //
