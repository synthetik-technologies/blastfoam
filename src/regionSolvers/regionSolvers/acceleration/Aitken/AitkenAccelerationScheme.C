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

#include "AitkenAccelerationScheme.H"
#include "valuePointPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::Aitken<Type, Patch, Mesh>::Aitken
(
    GeometricField<Type, Patch, Mesh>& field,
    autoPtr<PatchFieldSelector<Type>> selector,
    const dictionary& dict
)
:
    ResidualBase<Type, Patch, Mesh>(typeName, field, selector, dict),

    nFixed_(1),
    initRelaxFactor_(1.0),
    maxRelaxFactor_(1.0),
    aitkenFactor_(initRelaxFactor_)
{
    read(dict);
    Info<< indent << "nFixed: " << nFixed_ << nl
        << indent << "initialRelaxationFactor: " << initRelaxFactor_ << nl
        << indent << "maxRelaxationFactor: " << maxRelaxFactor_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::Aitken<Type, Patch, Mesh>::~Aitken()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::Aitken<Type, Patch, Mesh>::relax
(
    const label iter
)
{
    this->updateError();

    if (iter < 2)
    {
        aitkenFactor_ = initRelaxFactor_;
    }
    else
    {
        Field<Type> deltaResidual(this->residuals_ - this->prevResiduals_);
        scalar numerator =
            this->sumDotDot(this->prevResiduals_, deltaResidual);
        scalar denominator = this->sumDotDot(deltaResidual);

        if (denominator > small)
        {
            aitkenFactor_ = -aitkenFactor_*numerator/denominator;
        }
        else
        {
            aitkenFactor_ = maxRelaxFactor_;
        }

        if (mag(aitkenFactor_) > maxRelaxFactor_)
        {
            aitkenFactor_ = maxRelaxFactor_;
        }
    }

    Patch<Type>& pf = this->field_.boundaryFieldRef()[this->patchi_];
    Field<Type>& pfield =this->selector_->relaxField(pf);
    const Field<Type>& pfieldPrev =
        this->selector_->relaxField
        (
            this->field_.prevIter().boundaryField()[this->patchi_]
        );
    pfield = pfieldPrev + aitkenFactor_*this->residuals_;
    this->selector_->correct(pf);
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::Aitken<Type, Patch, Mesh>::read
(
    const dictionary& dict
)
{
    ResidualBase<Type, Patch, Mesh>::read(dict);
    this->coeffDict(dict).lookup("initialRelaxationFactor")
        >> initRelaxFactor_;
    this->coeffDict(dict).lookup("maxRelaxationFactor")
        >> maxRelaxFactor_;
    this->coeffDict(dict).readIfPresent("nFixed", nFixed_);
    if (nFixed_ < 1)
    {
        nFixed_ = max(nFixed_, 1);
    }
}


// ************************************************************************* //
