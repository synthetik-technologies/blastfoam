/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "SlaveLinearVelocityImmersedBoundaryObject.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::
SlaveLinearVelocityImmersedBoundaryObject
(
    const polyPatch& patch,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    ImmersedType(patch, dict, stateDict),
    patches_
    (
        dict.found("patches")
      ? dict.lookup("patches")
      : wordList(1, dict.lookup<word>("patch"))
    ),
    centerSet_(false),
    oldCentre_(Zero),
    centre_(Zero),
    velocity_(Zero)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::
~SlaveLinearVelocityImmersedBoundaryObject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
void
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::newTime()
{
    if (!this->pMesh_.moving())
    {
        velocity_ = Zero;
        return;
    }
    if (!centerSet_)
    {
        centre_ = this->shape().centre();
        oldCentre_ = centre_;
        movePoints();
    }

    const fvMesh& mesh = dynamicCast<const fvMesh>(this->pMesh_);
    const surfaceScalarField& meshPhi = mesh.phi();
    scalar maxp(max(this->points()).x());
    scalar maxfp;
    vector sumFlux(Zero);
    scalar sumArea(0.0);
    forAll(patches_, i)
    {
        const label patchID = mesh.boundaryMesh()[patches_[i]].index();
        const scalarField& magSf = mesh.magSf().boundaryField()[patchID];
        const vectorField& Sf = mesh.Sf().boundaryField()[patchID];
        maxfp = max(maxfp, max(mesh.C().boundaryField()[patchID]).x());
        sumFlux +=
            gSum
            (
                meshPhi.boundaryField()[patchID]*Sf/magSf
            );
        sumArea += gSum(magSf);
    }
    vector validD((vector(this->pMesh_.geometricD()) + vector::one)/2.0);
    velocity_ = vector(maxfp - maxp, 0, 0)/this->pMesh_.time().deltaTValue();//(sumFlux & vector(1, 0, 0))*vector(1, 0, 0)/sumArea ;

    oldCentre_ = centre_;
    centre_ += vector(maxfp - maxp, 0, 0);

    velocity_ = vector(centre_ - oldCentre_)/this->pMesh_.time().deltaTValue();
    Info<<centre_<<" "<<oldCentre_<<endl;;
Info<<maxp-maxfp<< " "<<velocity_<<endl;
    if (!centerSet_)
    {
        centerSet_ = true;
        velocity_ = Zero;
    }

    movePoints();
}


template<class ImmersedType>
void
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::update()
{
    ImmersedType::update();
}


template<class ImmersedType>
void
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::
movePoints()
{
    if (mag(centre_ - oldCentre_) < small)
    {
        return;
    }
    pointField newPoints(this->points() + (centre_ - oldCentre_));
    this->shape_->movePoints(newPoints);

    standAlonePatch::movePoints(newPoints);
    this->clearOut();
}


template<class ImmersedType>
void Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::status
(
    const bool print
) const
{
    ImmersedType::status(true);
    Info<< "    velocity: " << velocity_ << endl;
}


template<class ImmersedType>
Foam::point
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::centre() const
{
    return this->shape().centre();
}


template<class ImmersedType>
Foam::scalar
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::mass() const
{
    return great;
}


template<class ImmersedType>
Foam::diagTensor
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::momentOfInertia() const
{
    return diagTensor(great, great, great);
}


template<class ImmersedType>
Foam::vector
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::v
(
    const point& pt
) const
{
    return velocity_;
}


template<class ImmersedType>
Foam::tmp<Foam::vectorField>
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::velocity
(
    const pointField& points
) const
{
    return tmp<vectorField>(new vectorField(points.size(), velocity_));
}


template<class ImmersedType>
Foam::tmp<Foam::vectorField>
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::velocity() const
{
    return tmp<vectorField>(new vectorField(this->size(), velocity_));
}


template<class ImmersedType>
Foam::tensor
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::orientation() const
{
    return this->shape_->orientation();
}


template<class ImmersedType>
Foam::point
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::transform
(
    const point& initialPoint
) const
{
    return centre() + (orientation() & initialPoint);
}


template<class ImmersedType>
Foam::tmp<Foam::pointField>
Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::transform
(
    const pointField& initialPoints
) const
{
    return centre() + (orientation() & initialPoints);
}


template<class ImmersedType>
Foam::point Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::inverseTransform
(
    const point& pt
) const
{
    return (orientation().T() & (pt - centre()));
}


template<class ImmersedType>
Foam::tmp<Foam::pointField> Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::inverseTransform
(
    const pointField& points
) const
{
    return (orientation().T() & (points - centre()));
}


template<class ImmersedType>
bool Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::read
(
    const dictionary& dict
)
{
    return ImmersedType::read(dict);
}


template<class ImmersedType>
void Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::write
(
    Ostream& os
) const
{
    ImmersedType::write(os);
}


template<class ImmersedType>
void Foam::SlaveLinearVelocityImmersedBoundaryObject<ImmersedType>::write
(
    dictionary& dict
) const
{
    ImmersedType::write(dict);
}

// ************************************************************************* //
