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
    pBurst_(dict.lookup<scalar>("pBurst")),
    average_
    (
        !partialBurst_
      ? dict.lookup<bool>("useAverage")
      : dict.lookupOrDefault<bool>("useAverage", false)
    )
{
    if (average_ && partialBurst_)
    {
        WarningInFunction
            << "If partial burst is used, \"useAverage\" is ignored" << endl;
        average_ = false;
    }
}


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
        scalarField deltaP
        (
            this->patchField(p.boundaryField()[patch.index()])
        );
        if (!p.boundaryField()[patch.index()].coupled())
        {
            deltaP -= pRef_;
        }

        if (partialBurst_)
        {
            forAll(intact, facei)
            {
                if (deltaP[facei] > pBurst_)
                {
                    intact[facei] = 0;
                    burst = true;
                }
            }
            burst_ = gMax(intact) < small;
        }
        else if (average_)
        {
            burst =
                gSum(deltaP*patch.magSf())/gSum(patch.magSf()) > pBurst_;
            intact = !burst;
            burst_ = burst;
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
            << "Could not find " << pName_ << " field"
            << ", neglecting pressure" << endl
            << abort(FatalError);
    }
    return returnReduce(burst, orOp<bool>());
}


Foam::label Foam::burstModels::pressure::update
(
    const objectRegistry& obr,
    const labelList& ownCells,
    const labelList& neiCells,
    const scalarField& W
) const
{
    const volScalarField& p =
        obr.lookupObject<volScalarField>(pName_);
    scalar ownp = 0.0;
    scalar neip = 0.0;
    scalar ownW = 0.0;
    scalar neiW = 0.0;
    forAll(ownCells, i)
    {
        const label celli = ownCells[i];
        if (celli >= 0)
        {
            ownp += p[celli]*W[celli];
            ownW += W[celli];
        }
    }
    forAll(neiCells, i)
    {
        const label celli = neiCells[i];
        if (celli >= 0)
        {
            neip += p[celli]*W[celli];
            neiW += W[celli];
        }
    }
    reduce(ownp, sumOp<scalar>());
    reduce(ownW, sumOp<scalar>());
    reduce(neip, sumOp<scalar>());
    reduce(neiW, sumOp<scalar>());

    if (ownW > 0)
    {
        ownp /= ownW;
    }
    if (neiW > 0)
    {
        neip /= neiW;
    }

    if (mag(ownp - neip) > pBurst_)
    {
        return ownp > neip ? 1 : -1;
    }
    return false;
}


void Foam::burstModels::pressure::writeData(Ostream& os) const
{
    burstModel::writeData(os);
    writeEntry(os, "pName", pName_);
    writeEntry(os, "pBurst", pBurst_);
    writeEntry(os, "pRef", pRef_);
    writeEntry(os, "useAverage", average_);
}


// ************************************************************************* //
