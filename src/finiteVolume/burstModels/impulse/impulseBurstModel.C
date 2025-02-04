/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "impulseBurstModel.H"
#include "burstFvPatchFieldBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace burstModels
{
    defineTypeNameAndDebug(impulse, 0);
    addToRunTimeSelectionTable(burstModel, impulse, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModels::impulse::impulse(const dictionary& dict)
:
    burstModel(dict),
    impulseName_(dict.lookupOrDefault<word>("impulseName", "impulse")),
    impulseBurst_(dict.lookup<scalar>("impulseBurst"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModels::impulse::~impulse()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

bool Foam::burstModels::impulse::update
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
                if (deltaImp[facei] > impulseBurst_)
                {
                    intact[facei] = 0;
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
        WarningInFunction
            << "Could not find " << impulseName_
            << ", neglecting impulse. " << endl;
    }
    return returnReduce(burst, orOp<bool>());
}


Foam::label Foam::burstModels::impulse::update
(
    const objectRegistry& obr,
    const labelList& ownCells,
    const labelList& neiCells,
    const scalarField& W
) const
{
    const volScalarField& imp =
        obr.lookupObject<volScalarField>(impulseName_);
    scalar ownImp = 0.0;
    scalar neiImp = 0.0;
    scalar ownW = 0.0;
    scalar neiW = 0.0;
    forAll(ownCells, i)
    {
        const label celli = ownCells[i];
        if (celli >= 0)
        {
            ownImp += imp[celli]*W[celli];
            ownW += W[celli];
        }
    }
    forAll(neiCells, i)
    {
        const label celli = neiCells[i];
        if (celli >= 0)
        {
            neiImp += imp[celli]*W[celli];
            neiW += W[celli];
        }
    }
    reduce(ownImp, sumOp<scalar>());
    reduce(ownW, sumOp<scalar>());
    reduce(neiImp, sumOp<scalar>());
    reduce(neiW, sumOp<scalar>());

    if (ownW > 0)
    {
        ownImp /= ownW;
    }
    if (neiW > 0)
    {
        neiImp /= neiW;
    }

    if (log_)
    {
        Info<< indent << "Impulse differential: "<< mag(ownImp - neiImp)
            << endl;
    }

    if (mag(ownImp - neiImp) > impulseBurst_)
    {
        return ownImp > neiImp ? 1 : -1;
    }
    return false;
}


void Foam::burstModels::impulse::writeData(Ostream& os) const
{
    burstModel::writeData(os);
    writeEntry(os, "impulseName", impulseName_);
    writeEntry(os, "impulseBurst", impulseBurst_);
}


// ************************************************************************* //
