/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "MovingImmersedBoundaryObject.H"
#include "immersedBoundaryObjectMotionSolver.H"
#include "uniformDimensionedFields.H"
#include "septernion.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::applyRestraints()
{
    if (restraints_.empty())
    {
        return;
    }

    if (Pstream::master())
    {
        forAll(restraints_, rI)
        {
            if (report_)
            {
                Info<< "Restraint " << restraints_[rI].name() << ": ";
            }

            // Restraint position
            point rP = Zero;

            // Restraint force
            vector rF = Zero;

            // Restraint moment
            vector rM = Zero;

            // Accumulate the restraints
            restraints_[rI].restrain(*this, rP, rF, rM);

            // Update the acceleration
            a() += rF/this->mass_;

            // Moments are returned in global axes, transforming to
            // body local to add to torque.
            tau() += Q().T() & (rM + ((rP - centreOfRotation()) ^ rF));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::MovingImmersedBoundaryObject<ImmersedType>::
MovingImmersedBoundaryObject
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    ImmersedType(mesh, dict),
    motionState_(dict),
    motionState0_(),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialQ_
    (
        dict.lookupOrDefault
        (
            "initialOrientation",
            dict.lookupOrDefault("orientation", tensor::I)
        )
    ),
    momentOfInertia_(),
    aRelax_(dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0)),
    solver_(),
    points0_(this->points()),
    velocityBoundaries_(0)
{
    this->mass_ = readScalar(dict.lookup("mass"));

    // Scale the mass by the scale factor (axisymmetic)
    this->mass_ *= this->shape_->scale();
    momentOfInertia_ = this->mass()*this->shape_->momentOfInertia();


    this->initialCentreOfMass_ = dict.lookupOrDefault
    (
        "initialCentreOfMass",
        vector::zero
    );
    this->initialCentreOfRotation_ = dict.lookupOrDefault
    (
        "initialCentreOfRotation",
        this->initialCentreOfMass_
    );

    Q() &= this->shape_->orientation();

    centreOfRotation() += this->shape_->offset() + this->shape_->centreOfMass();
    centreOfMass() = centreOfRotation();
    solver_.set
    (
        immersedBoundaryObjectMotionSolver::New
        (
            dict.subDict("solver"),
            *this,
            motionState_,
            motionState0_
        ).ptr()
    );

    addRestraints(dict);

    // Set constraints and initial centre of rotation
    // if different to the centre of mass
    addConstraints(dict);

    // If the centres of mass and rotation are different ...
    vector R(this->initialCentreOfMass_ - this->initialCentreOfRotation_);
    if (magSqr(R) > vSmall)
    {
        // ... correct the moment of inertia tensor using parallel axes theorem
        momentOfInertia_ += this->mass_*diag(tensor::I*magSqr(R) - sqr(R));

        // ... and if the centre of rotation is not specified for motion state
        // update it
        if (!dict.found("centreOfRotation"))
        {
            motionState_.centreOfRotation() = this->initialCentreOfRotation_;
        }
    }

    // Save the old-time motion state
    motionState0_ = motionState_;

    movePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::MovingImmersedBoundaryObject<ImmersedType>::
~MovingImmersedBoundaryObject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::addRestraints
(
    const dictionary& dict
)
{
    if (dict.found("motionRestraints"))
    {
        const dictionary& restraintDict = dict.subDict("motionRestraints");

        label i = 0;

        restraints_.setSize(restraintDict.size());

        forAllConstIter(IDLList<entry>, restraintDict, iter)
        {
            if (iter().isDict())
            {
                restraints_.set
                (
                    i++,
                    immersedBoundaryObjectRestraint::New
                    (
                        this->pMesh_,
                        iter().keyword(),
                        iter().dict()
                    )
                );
            }
        }

        restraints_.setSize(i);
    }
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::addConstraints
(
    const dictionary& dict
)
{
    if (dict.found("motionConstraints"))
    {
        const dictionary& constraintDict = dict.subDict("motionConstraints");

        label i = 0;

        constraints_.setSize(constraintDict.size());

        pointConstraint pct;
        pointConstraint pcr;

        forAllConstIter(IDLList<entry>, constraintDict, iter)
        {
            if (iter().isDict())
            {
                constraints_.set
                (
                    i,
                    immersedBoundaryObjectConstraint::New
                    (
                        iter().keyword(),
                        iter().dict(),
                        *this
                    )
                );

                constraints_[i].setCentreOfRotation
                (
                    this->initialCentreOfRotation_
                );
                constraints_[i].constrainTranslation(pct);
                constraints_[i].constrainRotation(pcr);

                i++;
            }
        }

        constraints_.setSize(i);

        tConstraints_ = pct.constraintTransformation();
        rConstraints_ = pcr.constraintTransformation();

        Info<< "Translational constraint tensor " << tConstraints_ << nl
            << "Rotational constraint tensor " << rConstraints_ << endl;
    }
}


template<class ImmersedType>
Foam::Tuple2<Foam::tensor, Foam::vector>
Foam::MovingImmersedBoundaryObject<ImmersedType>::rotate
(
    const tensor& Q0,
    const vector& pi0,
    const scalar deltaT
) const
{
    Tuple2<tensor, vector> Qpi(Q0, pi0);
    tensor& Q = Qpi.first();
    vector& pi = Qpi.second();

    tensor R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorZ(deltaT*pi.z()/momentOfInertia_.zz());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    return Qpi;
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::updateAcceleration
(
    const vector& fGlobal,
    const vector& tauGlobal
)
{
    static bool first = true;

    // Save the previous iteration accelerations for relaxation
    vector aPrevIter = a();
    vector tauPrevIter = tau();

    // Calculate new accelerations
    a() = fGlobal/this->mass_;
    tau() = (Q().T() & tauGlobal);
    applyRestraints();

    // Relax accelerations on all but first iteration
    if (!first)
    {
        a() = aRelax_*a() + (1 - aRelax_)*aPrevIter;
        tau() = aRelax_*tau() + (1 - aRelax_)*tauPrevIter;
    }
    else
    {
        first = false;
    }
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::update()
{
    vectorField Md(this->faceCentres() - centreOfRotation());
    vectorField normal(this->Sf()/this->magSf());
    vectorField fN((this->force() & normal)*normal);
    vectorField fT(this->force() - fN);

    Pstream::listCombineGather(Md, plusEqOp<vector>());
    Pstream::listCombineGather(fN, plusEqOp<vector>());
    Pstream::listCombineGather(fT, plusEqOp<vector>());
    Pstream::listCombineScatter(Md);
    Pstream::listCombineScatter(fN);
    Pstream::listCombineScatter(fT);

    this->forceEff() += sum(fN) + sum(fT) + this->mass()*this->g_;
    this->momentEff() +=
        sum(Md ^ fN) + sum(Md ^ fT)
      + this->mass()*(momentArm() ^ this->g_);
    ImmersedType::update();

}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::solve(const bool restart)
{
    bool firstIter = true;

    if (restart)
    {
        motionState_ = motionState0_;
    }
    update
    (
        firstIter,
        this->forceEff() + this->forceExt(),
        this->momentEff() + this->momentExt(),
        this->pMesh_.time().deltaTValue(),
        this->pMesh_.time().deltaT0Value()
    );

    ImmersedType::solve();

    if (Pstream::master())
    {
        if (report_)
        {
            status();
        }
    }

}


template<class ImmersedType>
bool Foam::MovingImmersedBoundaryObject<ImmersedType>::update
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0
)
{
    if (Pstream::master())
    {
        solver_->solve(firstIter, fGlobal, tauGlobal, deltaT, deltaT0);
    }

    Pstream::scatter(motionState_);

    return true;
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::addVelocityBoundary
(
    volVectorField& U,
    const dictionary& defaultDict
)
{
    dictionary dict;
    bool found = false;
    if (this->dict_.found("boundaries"))
    {
        if (this->dict_.subDict("boundaries").found(U.name()))
        {
            dict = this->dict_.subDict("boundaries").subDict(U.name());
            found = true;
        }
    }
    if (!found)
    {
        dict = defaultDict.subDict(U.name());
    }
    velocityBoundaries_.resize(velocityBoundaries_.size() + 1);
    velocityBoundaries_.set
    (
        0,
        U.name(),
        immersedBoundaryVectorPatchField::New(U, dict, *this).ptr()
    );
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::addMomentumForcing
(
    const word& name,
    volVectorField& F,
    const volScalarField& alphaRho,
    const volVectorField& alphaRhoUOld,
    const volVectorField& RHS,
    const dimensionedScalar& dt
) const
{
    velocityBoundaries_[name].addForcing
    (
        F,
        alphaRho,
        alphaRhoUOld,
        RHS,
        dt.value()
    );
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::setValues()
{
    ImmersedType::setValues();
    forAll(velocityBoundaries_, i)
    {
        velocityBoundaries_[i].setValues();
    }
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::status() const
{
    ImmersedType::status();

    Info<< "    Centre of rotation: " << centreOfRotation() << nl
        << "    Centre of mass: " << centreOfMass() << nl
        << "    Orientation: " << orientation() << nl
        << "    Linear velocity: " << v() << nl
        << "    Angular velocity: " << omega()
        << endl;
}


template<class ImmersedType>
Foam::point Foam::MovingImmersedBoundaryObject<ImmersedType>::centreOfRotation() const
{
    return motionState_.centreOfRotation();
}


template<class ImmersedType>
Foam::point
Foam::MovingImmersedBoundaryObject<ImmersedType>::centreOfMass() const
{
    return transform(this->initialCentreOfMass_);
}


template<class ImmersedType>
Foam::tensor
Foam::MovingImmersedBoundaryObject<ImmersedType>::orientation() const
{
    return Q();
}


template<class ImmersedType>
Foam::vector
Foam::MovingImmersedBoundaryObject<ImmersedType>::momentArm() const
{
    return centreOfMass() - motionState_.centreOfRotation();
}


template<class ImmersedType>
Foam::vector
Foam::MovingImmersedBoundaryObject<ImmersedType>::v() const
{
    return motionState_.v();
}


template<class ImmersedType>
Foam::vector
Foam::MovingImmersedBoundaryObject<ImmersedType>::velocity
(
    const point& pt
) const
{
    return (omega() ^ (pt - centreOfRotation())) + v();
}


template<class ImmersedType>
Foam::tmp<Foam::vectorField>
Foam::MovingImmersedBoundaryObject<ImmersedType>::velocity
(
    const pointField& points
) const
{
    return (omega() ^ (points - centreOfRotation())) + v();
}


template<class ImmersedType>
Foam::vector
Foam::MovingImmersedBoundaryObject<ImmersedType>::omega() const
{
    return  Q() & (inv(momentOfInertia_) & pi());
}


template<class ImmersedType>
Foam::diagTensor
Foam::MovingImmersedBoundaryObject<ImmersedType>::momentOfInertia() const
{
    return momentOfInertia_;
}



template<class ImmersedType>
Foam::point
Foam::MovingImmersedBoundaryObject<ImmersedType>::transform
(
    const point& initialPoint
) const
{
    return
    (
        centreOfRotation()
      + (
            Q()
          & (initialPoint - this->initialCentreOfRotation_)
        )
    );
}


template<class ImmersedType>
Foam::tmp<Foam::pointField>
Foam::MovingImmersedBoundaryObject<ImmersedType>::transform
(
    const pointField& initialPoints
) const
{
    return
    (
        centreOfRotation()
      + (
            Q()
          & (initialPoints - this->initialCentreOfRotation_)
        )
    );
}


template<class ImmersedType>
Foam::tmp<Foam::pointField> Foam::MovingImmersedBoundaryObject<ImmersedType>::inverseTransform
(
    const pointField& points
) const
{
    tensor invQ(inv(Q()));
    return
    (
        (invQ & (points - centreOfRotation())) + this->initialCentreOfRotation_
    );
}


template<class ImmersedType>
Foam::tensor
Foam::MovingImmersedBoundaryObject<ImmersedType>::rConstraints() const
{
    return rConstraints_;
}


template<class ImmersedType>
Foam::tensor
Foam::MovingImmersedBoundaryObject<ImmersedType>::tConstraints() const
{
    return tConstraints_;
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::movePoints()
{
    ImmersedType::movePoints();
    this->shape_->movePoints();
    this->clearOut();
}


template<class ImmersedType>
bool Foam::MovingImmersedBoundaryObject<ImmersedType>::read
(
    const dictionary& dict
)
{
    this->mass_ = readScalar(dict.lookup("mass"));

    // Scale the mass by the scale factor (axisymmetic)
    this->mass_ *= this->shape_->scale();

    this->initialCentreOfMass_ = dict.lookupOrDefault
    (
        "initialCentreOfMass",
        vector(dict.lookup("centreOfMass"))
    );
    this->initialCentreOfRotation_ = this->initialCentreOfMass_;
    dict.lookup("momentOfInertia") >> momentOfInertia_;
    aRelax_ = dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0);

    restraints_.clear();
    addRestraints(dict);

    constraints_.clear();
    addConstraints(dict);

    return true;
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::write
(
    Ostream& os
) const
{
    motionState_.write(os);

    writeEntry(os, "mass", this->mass_/this->shape_->scale());
    writeEntry(os, "centreOfMass", this->initialCentreOfMass_);
    writeEntry(os, "initialOrientation", initialQ_);
    writeEntry(os, "momentOfInertia", momentOfInertia_);
    writeEntry(os, "accelerationRelaxation", aRelax_);
    writeEntry(os, "report", report_);

    if (!restraints_.empty())
    {
        os  << indent << "motionRestraints" << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;

        forAll(restraints_, rI)
        {
            word restraintType = restraints_[rI].type();

            os  << indent << restraints_[rI].name() << nl
                << indent << token::BEGIN_BLOCK << incrIndent << endl;

            writeEntry(os, "type", restraintType);

            restraints_[rI].write(os);

            os  << decrIndent << indent << token::END_BLOCK << endl;
        }

        os  << decrIndent << indent << token::END_BLOCK << nl;
    }

    if (!constraints_.empty())
    {
        os  << indent << "motionConstraints" << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;

        forAll(constraints_, rI)
        {
            word constraintType = constraints_[rI].type();

            os  << indent << constraints_[rI].name() << nl
                << indent << token::BEGIN_BLOCK << incrIndent << endl;

            writeEntry(os, "type", constraintType);

            constraints_[rI].immersedBoundaryObjectConstraint::write(os);

            constraints_[rI].write(os);

            os  << decrIndent << indent << token::END_BLOCK << endl;
        }

        os  << decrIndent << indent << token::END_BLOCK << nl;
    }
}


template<class ImmersedType>
void Foam::MovingImmersedBoundaryObject<ImmersedType>::write
(
    dictionary& dict
) const
{
    motionState_.write(dict);
    ImmersedType::write(dict);
}

// ************************************************************************* //
