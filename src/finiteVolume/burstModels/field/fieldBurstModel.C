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
#include "burstFvPatchFieldBase.H"
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

Foam::burstModels::field::field(const dictionary& dict)
:
    burstModel(dict),
    burstValues_()
{
    if (dict.found("values"))
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
    else
    {
        burstValues_ = HashTable<scalar>(dict.lookup("burstValues"));

    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burstModels::field::~field()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

bool Foam::burstModels::field::update
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
    forAllConstIter(HashTable<scalar>, burstValues_, iter)
    {
        const word& fieldName(iter.key());
        const scalar val = iter();
        if (mesh.foundObject<volScalarField>(fieldName))
        {
            const volScalarField& f =
                mesh.lookupObject<volScalarField>(fieldName);
            const fvPatchField<scalar>& pf =
                f.boundaryField()[patch.index()];
            scalarField deltaf(this->patchField(pf));

            if (partialBurst_)
            {
                forAll(intact, facei)
                {
                    if (deltaf[facei] > val)
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
                burst = gMax(deltaf) > val;
                intact = !burst;
                burst_ = burst;
            }
        }
        else
        {
            WarningInFunction
                << "Could not find " << fieldName
                << ", neglecting. " << endl;
        }
    }
    return returnReduce(burst, orOp<bool>());
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
