/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2022
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

#include "cellSizeErrorEstimator.H"
#include "meshSizeObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorEstimators
{
    defineTypeNameAndDebug(cellSize, 0);
    addToRunTimeSelectionTable(errorEstimator, cellSize, dictionary);
}
}

template<>
const char* Foam::NamedEnum<Foam::errorEstimators::cellSize::SizeType, 4>::names[] =
{
    "volume",
    "cmpt",
    "characteristic",
    "mag"
};

const Foam::NamedEnum<Foam::errorEstimators::cellSize::SizeType, 4>
Foam::errorEstimators::cellSize::sizeTypeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimators::cellSize::cellSize
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    errorEstimator(mesh, dict, name),
    sizeType_(CHARACTERISTIC),
    cmpts_()
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimators::cellSize::~cellSize()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::errorEstimators::cellSize::readCmpts(Istream& is) const
{
    labelHashSet cmpts;
    while (is.good())
    {
        token t(is);

        if (t.isNumber())
        {
            cmpts.insert(t.labelToken());
        }
        else if (t.isWord())
        {
            word c(t.wordToken());
            if (c == "x" || c == "X")
            {
                cmpts.insert(0);
            }
            else if (c == "y" || c == "Y")
            {
                cmpts.insert(1);
            }
            else if (c == "z" || c == "Z")
            {
                cmpts.insert(2);
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << "Unknown component " << c << endl
                    << "Valid components are: " << nl
                    << "x/X" << nl
                    << "y/Y" << nl
                    << "z/Z" << endl
                    << abort(FatalIOError);
            }
        }
    }
    const Vector<label> geoD(mesh_.geometricD());
    labelList cmptLables;
    forAllConstIter(labelHashSet, cmpts, iter)
    {
        if (geoD[iter.key()] == 1)
        {
            cmptLables.append(iter.key());
        }
    }
    return cmptLables;
}


Foam::labelList Foam::errorEstimators::cellSize::maxRefinement() const
{
    return mesh_.lookupObject<labelIOList>("cellLevel") + 1;
}


void Foam::errorEstimators::cellSize::update(const bool scale)
{
    if (updateCurTimeIndex(!scale))
    {
        return;
    }
    const labelHashSet eCells(this->errorCells());

    scalarField& errorI = error_;
    if (sizeType_ == VOLUME)
    {
        const scalarField& V = mesh_.V();
        forAllConstIter(labelHashSet, eCells, iter)
        {
            const label celli = iter.key();
            errorI[celli] = V[celli];
        }
    }
    else if (sizeType_ == CHARACTERISTIC)
    {
        const scalarField& dx = meshSizeObject::New(mesh_).dx();
        forAllConstIter(labelHashSet, eCells, iter)
        {
            const label celli = iter.key();
            errorI[celli] = dx[celli];
        }
    }
    else if (sizeType_ == MAG)
    {
        const vectorField& dX = meshSizeObject::New(mesh_).dX();
        const Vector<label>& geoD = mesh_.geometricD();
        forAllConstIter(labelHashSet, eCells, iter)
        {
            const label celli = iter.key();
            forAll(geoD, cmpti)
            {
                if (geoD[cmpti] > 0)
                {
                    errorI[celli] += sqr(component(dX[celli], cmpti));
                }
            }
            errorI[celli] = sqrt(errorI[celli]);
        }
    }
    else if (sizeType_ == CMPT)
    {
        const vectorField& dX = meshSizeObject::New(mesh_).dX();
        forAllConstIter(labelHashSet, eCells, iter)
        {
            const label celli = iter.key();
            forAll(cmpts_, i)
            {
                const label cmpti = cmpts_[i];
                if (dX[celli][cmpti] > maxDX_[cmpti])
                {
                    errorI[celli] = max(1.0, errorI[celli]);
                }
                else if (dX[celli][cmpti] < minDX_[cmpti])
                {
                    errorI[celli] = max(-1.0, errorI[celli]);
                }
                else
                {
                    errorI[celli] = max(0, errorI[celli]);
                }
            }
        }
    }

    if (sizeType_ != CMPT)
    {
        forAllConstIter(labelHashSet, eCells, iter)
        {
            const label celli = iter.key();
            if (errorI[celli] < lowerUnrefine_)
            {
                errorI[celli] = -1.0;
            }
            else if (errorI[celli] > lowerRefine_)
            {
                errorI[celli] = 1.0;
            }
            else
            {
                errorI[celli] = 0.0;
            }
        }
    }
}


void Foam::errorEstimators::cellSize::read(const dictionary& dict)
{
    sizeType_ = sizeTypeNames.read(dict.lookup("type"));
    if (sizeType_ == CMPT)
    {
        cmpts_ = readCmpts(dict.lookup("cmpts"));
        minDX_ = dict.lookup<vector>("minDX");
        maxDX_ = dict.lookup<vector>("maxDX");
    }
    else if (sizeType_ == VOLUME)
    {
        lowerUnrefine_ = dict.lookup<scalar>("minVolume");
        lowerRefine_ = dict.lookup<scalar>("maxVolume");
    }
    else
    {
        lowerUnrefine_ = dict.lookup<scalar>("minDx");
        lowerRefine_ = dict.lookup<scalar>("maxDx");
    }

    maxLevel_ = dict.lookupOrDefault<label>("maxRefinement", 10);

    readCellZones(dict);
}

// ************************************************************************* //
