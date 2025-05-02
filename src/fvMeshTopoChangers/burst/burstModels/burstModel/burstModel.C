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

#include "burstModel.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstModel, 0);
    defineRunTimeSelectionTable(burstModel, dictionary);
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::burstModel> Foam::burstModel::New
(
    const dictionary& dict,
    const bool coupled
)
{
    word modelType(dict.lookup("burstModel"));

    Info<< "Selecting burst model: " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown burst model "
            << modelType << endl << endl
            << "Valid burst models are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, coupled);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModel::burstModel
(
    const dictionary& dict,
    const bool coupled
)
:
    partialBurst_(dict.lookup<bool>("partialBurst")),
    useDelta_
    (
        coupled
      ? dict.lookupOrDefault<bool>("useDelta", true)
      : false
    ),
    useAverage_(dict.lookupOrDefault<bool>("useAverage", false)),
    log_(dict.lookupOrDefault("logBurst", false)),
    needUpdate_(true)
{

    if (useAverage_ && partialBurst_)
    {
        WarningInFunction
            << "If partial burst is used, \"useAverage\" is ignored" << endl;
        useAverage_ = false;
    }

    if (coupled)
    {
        mappingPtr_ =
            patchToPatch::New
            (
                dict.lookup<word>("patchToPatch"),
                false
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModel::~burstModel()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::burstModel::needUpdate() const
{
    needUpdate_ = true;
}


void Foam::burstModel::updateMapping
(
    const polyPatch& patch1,
    const polyPatch& patch2
) const
{
    if (!needUpdate_)
    {
        return;
    }
    mappingPtr_->update
    (
        patch1,
        PatchTools::pointNormals(patch1.boundaryMesh().mesh(), patch1),
        patch2
    );
    needUpdate_ = false;
}


bool Foam::burstModel::findBurstFaces
(
    const scalarField& pf1,
    const fvPatch&  patch1,
    const scalarField& pf2,
    const fvPatch&  patch2,
    const scalar burstValue,
    labelHashSet& master,
    labelHashSet& slave
) const
{
    bool allBurst = false;
    const polyPatch& pp1 = patch1.patch();
    const polyPatch& pp2 = patch2.patch();

    updateMapping(pp1, pp2);

    if (useAverage_)
    {
        const scalar pf1Mean = gSum(pf1*patch1.magSf())/gSum(patch1.magSf());
        const scalar pf2Mean = gSum(pf2*patch2.magSf())/gSum(patch2.magSf());
        allBurst =
            (useDelta_ && mag(pf1Mean - pf2Mean) > burstValue)
         || (!useDelta_ && max(pf1Mean, pf2Mean) > burstValue);
    }
    else
    {
        // Map pf1 and pf2 to the other side
        scalarField refVal1
        (
            useDelta_
          ? mag(pf1 - mappingPtr_->tgtToSrc(pf2, pf1))
          : max(pf1,  mappingPtr_->tgtToSrc(pf2, pf1))
        );
        refVal1 -= burstValue;

        scalarField refVal2
        (
            useDelta_
          ? mag(pf2 - mappingPtr_->tgtToSrc(pf1, pf2))
          : max(pf2,  mappingPtr_->tgtToSrc(pf1, pf2))
        );
        refVal2 -= burstValue;

        if (partialBurst_)
        {
            forAll(pp1, fi)
            {
                if (refVal1[fi] > 0)
                {
                    master.insert(pp1.start() + fi);
                }
            }

            forAll(pp2, fi)
            {
                if (refVal2[fi] > 0)
                {
                    slave.insert(pp2.start() + fi);
                }
            }
        }
        else
        {
            allBurst = gMax(refVal1) > 0 || gMax(refVal2) > 0;
        }
    }
    reduce(allBurst, orOp<bool>());


    if (allBurst)
    {
        forAll(pp1, fi)
        {
            master.insert(pp1.start() + fi);
        }

        forAll(pp2, fi)
        {
            slave.insert(pp2.start() + fi);
        }
    }
    return allBurst;
}


void Foam::burstModel::writeData(Ostream& os) const
{
    writeEntry(os, "burstModel", type());
    writeEntry(os, "partialBurst", partialBurst_);
    writeEntry(os, "useDelta", useDelta_);
    writeEntry(os, "useAverage", useAverage_);
    writeEntryIfDifferent(os, "log", log_, false);
}


// ************************************************************************* //
