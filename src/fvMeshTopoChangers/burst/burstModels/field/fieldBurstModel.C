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

#include "fieldBurstModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace burstModels
{
    defineTypeNameAndDebug(field, 0);
    addToRunTimeSelectionTable(burstModel, field, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burstModels::field::field
(
    const dictionary& dict,
    const bool coupled
)
:
    burstModel(dict, coupled),
    burstValues_()
{
    if
    (
        dict.found("values")
     && dict.found("fields")
    )
    {
        List<word> fields(dict.lookup("fields"));
        List<scalar> values(dict.lookup("values"));
        forAll(fields, i)
        {
            burstValues_.insert(fields[i], values[i]);
        }
    }
    else if (dict.found("burstValue"))
    {
        Tuple2<word, scalar> nameVal(dict.lookup("burstValue"));
        burstValues_.insert(nameVal.first(), nameVal.second());
    }
    else if (dict.found("burstValues"))
    {
        burstValues_ = HashTable<scalar>(dict.lookup("burstValues"));
    }
    else
    {
        FatalErrorInFunction
            << "Either fields and values, burstValue, or burstValues "
            << "must be provided" << endl
            << abort(FatalError);
    }

    // if (log_)
    // {
    //     osPtr_() << '#' << '\t' << "time" << '\t';
    //     forAllConstIter(HashTable<scalar>, burstValues_, iter)
    //     {
    //         osPtr_() << iter.key() << '\t';
    //     }
    //     osPtr_() << endl;
    // }
}


Foam::burstModels::field::field
(
    const dictionary& dict,
    const bool coupled,
    const wordList& names,
    const scalarList& vals
)
:
    burstModel(dict, coupled),
    burstValues_()
{

    forAll(names, i)
    {
        burstValues_.insert(names[i], vals[i]);
    }
    // if (log_)
    // {
    //     osPtr_() << '#' << '\t' << "time" << '\t';
    //     forAllConstIter(HashTable<scalar>, burstValues_, iter)
    //     {
    //         osPtr_() << iter.key() << '\t';
    //     }
    //     osPtr_() << endl;
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModels::field::~field()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

bool Foam::burstModels::field::facesToChange
(
    const fvPatch& patch1,
    const fvPatch& patch2,
    labelList& masterFaces,
    labelList& slaveFaces
) const
{
    if
    (
        !returnReduce(patch1.size(), sumOp<label>())
     || !returnReduce(patch2.size(), sumOp<label>())
    )
    {
        return false;
    }

    const fvMesh& mesh = patch1.boundaryMesh().mesh();
    labelHashSet masterSet, slaveSet;

    forAllConstIter(HashTable<scalar>, burstValues_, iter)
    {
        const word& fieldName(iter.key());
        if (mesh.foundObject<volScalarField>(fieldName))
        {
            const volScalarField& f =
                mesh.lookupObject<volScalarField>(fieldName);
            if
            (
                this->findBurstFaces
                (
                    f.boundaryField()[patch1.index()],
                    patch1,
                    f.boundaryField()[patch2.index()],
                    patch2,
                    iter(),
                    masterSet,
                    slaveSet
                )
            )
            {
                break;
            }
        }
        else
        {
            WarningInFunction
                << "Could not find " << fieldName
                << ", neglecting. " << endl;
        }
    }
    masterFaces = masterSet.toc();
    slaveFaces = slaveSet.toc();
    return returnReduce(masterFaces.size() + slaveFaces.size(), orOp<bool>());
}


bool Foam::burstModels::field::facesToChange
(
    const fvPatch& patch,
    labelList& faces
) const
{
    if (!returnReduce(patch.size(), sumOp<label>()))
    {
        return false;
    }

    const fvMesh& mesh = patch.boundaryMesh().mesh();
    const polyPatch& pp = patch.patch();
    labelHashSet faceSet;

    forAllConstIter(HashTable<scalar>, burstValues_, iter)
    {
        const word& fieldName(iter.key());
        if (mesh.foundObject<volScalarField>(fieldName))
        {
            const volScalarField& f =
                mesh.lookupObject<volScalarField>(fieldName);
            const scalarField& pf = f.boundaryField()[patch.index()];

            bool allBurst = false;

            if (useAverage_)
            {
                const scalar pfMean =
                    gSum(pf*patch.magSf())/gSum(patch.magSf());
                allBurst = pfMean > iter();
            }
            else
            {
                // Map pf1 and pf2 to the other side
                const scalarField refVal(pf - iter());

                if (partialBurst_)
                {
                    forAll(pp, fi)
                    {
                        if (refVal[fi] > 0)
                        {
                            faceSet.insert(pp.start() + fi);
                        }
                    }
                }
                else
                {
                    allBurst = gMax(refVal) > 0;
                }
            }
            reduce(allBurst, orOp<bool>());


            if (allBurst)
            {
                forAll(pp, fi)
                {
                    faceSet.insert(pp.start() + fi);
                }
                break;
            }
        }
        else
        {
            WarningInFunction
                << "Could not find " << fieldName
                << ", neglecting. " << endl;
        }
    }
    faces = faceSet.toc();
    return returnReduce(faces.size(), orOp<bool>());
}


Foam::label Foam::burstModels::field::update
(
    const objectRegistry& obr,
    const labelList& ownCells,
    const labelList& neiCells,
    const scalarField& W
) const
{
    scalar ownW = 0.0;
    scalar neiW = 0.0;
    forAll(ownCells, i)
    {
        if (ownCells[i] >= 0)
        {
            ownW += W[ownCells[i]];
        }
    }
    forAll(neiCells, i)
    {
        if (neiCells[i] >= 0)
        {
            neiW += W[neiCells[i]];
        }
    }
    reduce(ownW, sumOp<scalar>());
    reduce(neiW, sumOp<scalar>());

    forAllConstIter(HashTable<scalar>, burstValues_, iter)
    {
        const word& fieldName(iter.key());
        const scalar val = iter();

        const volScalarField& f =
            obr.lookupObject<volScalarField>(fieldName);

        scalar ownF = 0.0;
        scalar neiF = 0.0;
        forAll(ownCells, i)
        {
            const label celli = ownCells[i];
            if (celli >= 0)
            {
                ownF += f[celli]*W[celli];
            }
        }
        forAll(neiCells, i)
        {
            const label celli = neiCells[i];
            if (celli >= 0)
            {
                neiF += f[celli]*W[celli];
            }
        }
        reduce(ownF, sumOp<scalar>());
        reduce(neiF, sumOp<scalar>());

        if (ownW > 0)
        {
            ownF /= ownW;
        }
        if (neiW > 0)
        {
            neiF /= neiW;
        }

        if (log_)
        {
            Info<< indent << iter.key() << " differential: "<< mag(ownF - neiF)
                << endl;
        }
        if (mag(ownF - neiF) > val)
        {
            return ownF > neiF ? 1 : -1;
        }
    }
    return false;
}


void Foam::burstModels::field::writeData(Ostream& os) const
{
    burstModel::writeData(os);

    // Write as two seperate lists to remove errors in Paraview
    List<word> fields(burstValues_.toc());
    List<scalar> vals(fields.size());
    forAll(fields, i)
    {
        vals[i] = burstValues_[fields[i]];
    }
    writeEntry(os, "fields", fields);
    writeEntry(os, "values", vals);
}


// ************************************************************************* //
