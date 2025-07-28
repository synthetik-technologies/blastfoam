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

#include "IQNILSAccelerationScheme.H"
#include "scalarMatrices.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::IQNILS<Type, Patch, Mesh>::GramSchmidt_QRdecomp
(
    Field<Type>& pfieldNew
)
{
    // Calculate the QR-decomposition of V
    // with modified Gram-Schmidt
    label n = this->Vs_.size();
    scalarSquareMatrix R(n, 0.0);
    Field<scalar> alpha(n, Zero);
    List<Field<Type>> Q(n);
    for (label i = 0; i < n; i++)
    {
        Q[i] = this->Vs_[n-1-i];
    }

    for (label i = 0; i < n; i++)
    {
        // Project proj_v(u) = <u, v>/<v, v>

        // Magnitude
        R[i][i] = sqrt(this->sumDotDot(Q[i]));

        // Make sure magnitude is not negligible
        // Ignore contribution otherwise (alpha = 0)
        if (R[i][i] > small)
        {
            Q[i] /= R[i][i];

            for (label j = i+1; j < n; j++)
            {
                // Dot product
                R[i][j] = this->sumDotDot(Q[i], Q[j]);
                Q[j] -= R[i][j]*Q[i];
            }

            alpha[i] = -this->sumDotDot(Q[i], this->residuals_);
        }
        else
        {
            R[i][i] = 0.0;
        }
    }

    scalar epsilon = 0.0;
    for (label j = 0; j < n; j++)
    {
        scalar RColSum = 0.0;
        for (label i = 0; i < j+1; i++)
        {
            RColSum += mag(R[i][j]);
        }
        epsilon = max(epsilon, RColSum);
    }
    epsilon *= 1e-10;

    for (label i = 0; i < n; i++)
    {
        if (mag(R[i][i]) > epsilon)
        {
            for (label j = i+1; j < n; j++)
            {
                R[i][j] /= R[i][i];
            }
            alpha[i] /= R[i][i];
            R[i][i] = 1.0;
        }
    }

    for (label j = n-1; j >= 0; j--)
    {
        if (mag(R[j][j]) > epsilon)
        {
            for (label i = 0; i < j; i++)
            {
                alpha[i] -= alpha[j]*R[i][j];
            }
        }
        else
        {
            alpha[j] = 0.0;
        }
    }

    forAll(this->Ws_, i)
    {
        if (mag(alpha[n-1-i]) > small)
        {
            // // Limit alpha
            // if (mag(alpha[n-1-i]) > 2.0)
            // {
            //     alpha[n-1-i] = sign(alpha[n-1-i]);
            // }
            pfieldNew += this->Ws_[i]*alpha[n-1-i];
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::IQNILS<Type, Patch, Mesh>::IQNILS
(
    GeometricField<Type, Patch, Mesh>& field,
    autoPtr<PatchFieldSelector<Type>> selector,
    const dictionary& dict
)
:
    QNBase<Type, Patch, Mesh>(typeName, field, selector, dict),

    initRelaxFactor_(1.0),
    nFixed_(1),
    qr_("GramSchmidt")
{
    read(dict);
    Info<< indent << "initialRelaxationFactor: " << initRelaxFactor_ << nl
        << indent << "nFixed: " << nFixed_ << nl
        << indent << "QRDeomposition: " << qr_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::IQNILS<Type, Patch, Mesh>::~IQNILS()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::IQNILS<Type, Patch, Mesh>::relax
(
    const label iter
)
{
    this->updateError();

    Patch<Type>& pf = this->field_.boundaryFieldRef()[this->patchi_];
    Field<Type>& pfield =this->selector_->relaxField(pf);
    const Field<Type>& pfieldPrev =
        this->selector_->relaxField
        (
            this->field_.prevIter().boundaryField()[this->patchi_]
        );
    const Field<Type> oldUnrelaxed(pfield);


    if (iter == nFixed_)
    {
        // Set the reference state
        this->updateVW(0);
    }
    else if (iter > nFixed_)
    {
        // Append residuals
        this->updateVW(1);
    }

    //- Set the reference field to the current unrelaxed field
    this->refField_ = oldUnrelaxed;

    if (iter < 0)
    {
        // No relaxation
    }
    else if (iter < nFixed_ || this->times_.size() < 2)
    {
        pfield = pfieldPrev + initRelaxFactor_*this->residuals_;
    }
    else
    {
        Field<Type> pfieldNew(dynamicCast<const Field<Type>>(pfield));

        if (qr_ == "GramSchmidt")
        {
            GramSchmidt_QRdecomp(pfieldNew);
        }
        else
        {
            FatalErrorInFunction
            << "Unknown QR decomposition method " << qr_ << nl
            << "valid options are:" << nl
            << "GramSchmidt" << nl
            << endl
            << abort(FatalError);
        }

        pfield = pfieldNew;
    }
    this->selector_->correct(pf);
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::IQNILS<Type, Patch, Mesh>::read
(
    const dictionary& dict
)
{
    QNBase<Type, Patch, Mesh>::read(dict);
    this->coeffDict(dict).lookup("initialRelaxationFactor")
        >> initRelaxFactor_;
    this->coeffDict(dict).readIfPresent("nFixed", nFixed_);
    if (nFixed_ < 1)
    {
        nFixed_ = max(nFixed_, 1);
    }

    // this->coeffDict(dict).readIfPresent("QRDecomposition", qr_);
    // if (qr_ != "GramSchmidt")
    // {
    //     FatalErrorInFunction
    //         << "Unknown QR decomposition method " << qr_ << nl
    //         << "valid options are:" << nl
    //         << "GramSchmidt" << nl
    //         << endl
    //         << abort(FatalError);
    // }
}


// ************************************************************************* //
