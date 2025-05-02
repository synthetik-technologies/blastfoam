/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "tractionBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tractionBase::tractionBase(const fvPatch& p)
:
    patch_(p),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    force_(Zero)
{}


Foam::tractionBase::tractionBase
(
    const fvPatch& p,
    const dictionary& dict,
    const bool mustRead
)
:
    patch_(p),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    force_(Zero)
{
    if (mustRead)
    {
        traction_ = vectorField("traction", dict, p.size());
        pressure_ = scalarField("pressure", dict, p.size());
    }
}


Foam::tractionBase::tractionBase
(
    const tractionBase& tb,
    const fvPatch& p,
    const fieldMapper& mapper
)
:
    patch_(p),
    traction_(p.size(), Zero),
    pressure_(p.size(), Zero),
    force_(tb.force_)
{
    mapper(traction_, tb.traction_);
    mapper(pressure_, tb.pressure_);
}


Foam::tractionBase::tractionBase
(
    const tractionBase& tb,
    const fvPatch& p
)
:
    patch_(p),
    traction_(tb.traction_),
    pressure_(tb.pressure_),
    force_(tb.force_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tractionBase::map
(
    const fvPatchField<vector>& ptf,
    const fieldMapper& mapper
)
{
    const tractionBase& tb = dynamicCast<const tractionBase>(ptf);
    mapper(traction_, tb.traction_);
    mapper(pressure_, tb.pressure_);
}


void Foam::tractionBase::reset(const fvPatchField<vector>& ptf)
{
    const tractionBase& tb = dynamicCast<const tractionBase>(ptf);
    traction_.reset(tb.traction_);
    pressure_.reset(tb.pressure_);
}


void Foam::tractionBase::updateForce()
{
    force_ =
        gSum
        (
            pressure_*patch_.Sf()
          + traction_*patch_.magSf()
        );
}


void Foam::tractionBase::write(Ostream& os) const
{
    writeEntry(os, "traction", traction_);
    writeEntry(os, "pressure", pressure_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tractionBase, 0);
}

bool Foam::tractionBase::canRelax = true;


// ************************************************************************* //
