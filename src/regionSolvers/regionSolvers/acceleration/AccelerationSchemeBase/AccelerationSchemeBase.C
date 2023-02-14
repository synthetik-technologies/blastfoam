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

#include "AccelerationSchemeBase.H"
#include "PatchFieldSelector.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::AccelerationSchemeBase<Type, Patch, Mesh>::AccelerationSchemeBase
(
    const word& type,
    GeometricField<Type, Patch, Mesh>& field,
    autoPtr<PatchFieldSelector<Type>> selector,
    const dictionary& dict
)
:
    accelerationScheme(type, selector->index(), dict),
    field_(field),
    selector_(selector)
{
    regionName_ = field.mesh().thisDb().name();
    fieldName_ = field.name();
    patchName_= field.mesh().boundary()[selector_->index()].name();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::AccelerationSchemeBase<Type, Patch, Mesh>::~AccelerationSchemeBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::AccelerationSchemeBase<Type, Patch, Mesh>::updateError()
{
    error_ = Zero;

    const Field<Type>& pfield =
        selector_->relaxField(field_.boundaryField()[patchi_]);
    const Field<Type>& pfieldPrev =
        selector_->relaxField(field_.prevIter().boundaryField()[patchi_]);

    scalar maxErrorMagSqr = 0.0;
    const label n = returnReduce(pfield.size(), sumOp<label>());

    forAll(pfield, i)
    {
        scalar errorMagSqr = magSqr(pfield[i] - pfieldPrev[i]);
        error_ += errorMagSqr;
        maxErrorMagSqr = max(maxErrorMagSqr, errorMagSqr);
    }

    reduce(error_, sumOp<scalar>());
    reduce(maxErrorMagSqr, maxOp<scalar>());
    error_ = sqrt(error_/(maxErrorMagSqr + small))/scalar(n);

    if (initialError_ < 0)
    {
        initialError_ = error_;
    }
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::AccelerationSchemeBase<Type, Patch, Mesh>::storePrevIter()
{
    const Patch<Type>& pfield = field_.boundaryField()[patchi_];
    Patch<Type>& pfieldPrev =
        const_cast<Patch<Type>&>(field_.prevIter().boundaryField()[patchi_]);
    selector_->storePrevIter(pfield, pfieldPrev);
}


// ************************************************************************* //
