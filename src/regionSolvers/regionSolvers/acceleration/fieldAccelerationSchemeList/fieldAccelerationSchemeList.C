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

#include "fieldAccelerationSchemeList.H"
#include "globalPolyBoundaryMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldAccelerationSchemeList, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldAccelerationSchemeList::fieldAccelerationSchemeList
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& subDictName
)
:
    PtrListDictionary<fieldAccelerationScheme>(0),
    mesh_(mesh),
    dict_(dict),
    subDictName_(subDictName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldAccelerationSchemeList::~fieldAccelerationSchemeList()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldAccelerationSchemeList::addField(const word& fieldName)
{
    addField(fieldName, globalPolyBoundaryMesh::New(mesh_).coupledPatches());
}


void Foam::fieldAccelerationSchemeList::addField
(
    const word& fieldName,
    const labelList& patches
)
{
    if (this->found(fieldName))
    {
        return;
    }

    label fieldi = this->size();
    this->setSize(fieldi + 1);
    this->set
    (
        fieldi,
        fieldName,
        new fieldAccelerationScheme
        (
            fieldName,
            mesh_,
            patches,
            (subDictName_ != word::null)
          ? dict_.subDict(subDictName_)
          : dict_
        )
    );
}


void Foam::fieldAccelerationSchemeList::updateError()
{
    forAll(*this, i)
    {
        this->operator[](i).updateError();
    }
}


void Foam::fieldAccelerationSchemeList::relax(const label iter)
{
    forAll(*this, i)
    {
        this->operator[](i).relax(iter);
    }
}


void Foam::fieldAccelerationSchemeList::clear(const bool full)
{
    forAll(*this, i)
    {
        this->operator[](i).clear(full);
    }
}


Foam::Convergence Foam::fieldAccelerationSchemeList::convergence() const
{
    bool allUnknown = true;
    bool allFull = true;
    forAll(*this, fieldi)
    {
        switch (operator[](fieldi).convergence())
        {
            case NOT_CONVERGED:
                return NOT_CONVERGED;
            case CONVERGED:
                allFull = false;
                allUnknown = false;
                break;
            case FULL_CONVERGENCE:
                allUnknown = false;
                break;
            default:
                break;
        }
    }
    return
        allUnknown
      ? UNKNOWN_CONVERGENCE
      : (
            allFull ? FULL_CONVERGENCE : CONVERGED
        );
}


void Foam::fieldAccelerationSchemeList::storePrevIter()
{
    forAll(*this, i)
    {
        this->operator[](i).storePrevIter();
    }
}


bool Foam::fieldAccelerationSchemeList::converged() const
{
    Convergence converged = this->convergence();
    return
        converged == CONVERGED
     || converged == FULL_CONVERGENCE;
}


void Foam::fieldAccelerationSchemeList::read(const dictionary& dict)
{
    forAll(*this, i)
    {
        this->operator[](i).read(dict);
    }
}


void Foam::fieldAccelerationSchemeList::print(Ostream& os) const
{
    forAll(*this, i)
    {
        this->operator[](i).print(os);
    }
}


// ************************************************************************* //
