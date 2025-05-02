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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "IterativeSolidModel.H"

// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

void Foam::QuasiNewtonSolidRelaxationModel::readDict()
{
    const dictionary& dict = this->coeffDict();
    dict.readIfPresent("restartFrequency", restartFreq_);

    V_.setSize(restartFreq_ + 2),
    W_.setSize(restartFreq_ + 2),
    T_.setSize(rRestartFreq_ + 2),
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::QuasiNewtonSolidRelaxationModel::QuasiNewtonSolidRelaxationModel
(
    solidModel& solid,
)
:
    solidRelaxationModel(solid),
    restartFreq_(25),
    V_(),
    W_(),
    T_(),
    DRef_
    (
        IOobject
        (
            solid_.solutionD().name() + "Ref",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    unrelaxedDRef_
    (
        IOobject
        (
            "unrelaxted" + solid_.solutionD().name() + "Ref",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    unrelaxedDRef_()
{
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class SolidModel>
Foam::IterativeSolidModel<SolidModel>::~IterativeSolidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::QuasiNewtonSolidRelaxationModel::relaxField
(
    volVectorField& D,
    const label iCorr
)
{

    // This method is a modified form of the IQNILS by Degroote et al.

    // J. Degroote, K.-J. Bathe and J. Vierendeels.
    // A fluid solid interaction solver with IQN-ILS coupling algorithm.
    // Performance of a new partitioned procedure versus a monolithic
    // procedure in fluid-solid interaction. Computers & Solids

    if (iCorr == 0 || iCorr % restartFreq_ == 0)
    {
        // Clean up data from old time steps

        if (solidModel::debug)
        {
            Info<< "Modes before clean-up : " << T_.size();
        }

        while (true)
        {
            if (T_.size())
            {
                if
                (
                    this->runTime().timeIndex() > T_[0]
                 || (iCorr % restartFreq_) == 0
                )
                {
                    for (label i = 0; i < T_.size() - 1; i++)
                    {
                        T_[i] = T_[i + 1];
                        V_[i].transfer(V_[i + 1]);
                        W_[i].transfer(W_[i + 1]);
                    }

                    T_.remove();
                    V_.remove();
                    W_.remove();
                }
                else
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }

        if (solidModel::debug)
        {
            Info<< ", modes after clean-up : " << T_.size() << endl;
        }
    }
    else if (iCorr == 1 || iCorr % restartFreq_ == 1)
    {
        // Set reference in the first coupling iteration
        unrelaxedDRef_ == D;
        DRef_ == D.prevIter();
    }
    else
    {
        // Store the input vector field, defined as the previous iteration
        // D field (after relaxation) minus the Dp previous iteration field
        // in the first iteration (after relaxation)
        V_.append
        (
            (D - D.prevIter()) - (unrelaxedDRef_ - DRef_)
        );

        // V should be (from FSI paper):
        // DeltaR^{k-1} = R^{k-1} - R^k = DPrevIt.PrevIt - DPrevIter
        // DeltaR^{k-2} = R^{k-2} - R^k = D.PI.PI.PI - D.PI
        // ...
        // DeltaR^{0} =  R^0 - R^k = DRef - D.prevIter
        // Or in the general paper:
        // V_i = p_k - p_i    for i = 0, 1, ..., k - 1
        // V_{k-1} = p_k - p_{k-1} = D.PI - D.PI.PI
        // V_{k-2} = p_k - p_{k-2} = D.PI - D.PI.PI.PI
        // ...
        // V_{0} = p_k - p_{0} = D.PI - DRef
        // BUT, the implemented code does this:
        // V_i = p_{i+1} - p_0    for i = 0, 1, ..., k - 1
        // V_{k-1} = p_{k} - p_0 = D.PI - DRef
        // V_{k-2} = p_{k-1} - p_0 = D.PI{k-1} - DRef
        // ...
        // V_{0} = p_{1} - p_{0} = D.PI_1 - DRef
        // This means that we just append  the following line each
        // iteration:
        // V_{k-1} = p_{k} - p_0 = D.PI - DRef
        // It this equivalent?
        // We could try implementing it as described in the paper, but this
        // will require D and D.prevIter and their history

        // Store the output vector field, defined as the current iteration
        // D field (before relaxation) minus the D field in the first
        // iteration (before relaxation)
        W_.append(D - unrelaxedDRef_);

        // Store the time index
        T_.append(this->runTime().timeIndex());
    }

    if (T_.size() > 1)
    {
        // Consider QuasiNewtonV as a matrix V
        // with as columns the items
        // in the DynamicList and calculate the QR-decomposition of V
        // with modified Gram-Schmidt
        label cols = V_.size();
        RectangularMatrix<scalar> R(cols, cols, 0.0);
        RectangularMatrix<scalar> C(cols, 1);
        RectangularMatrix<scalar> Rcolsum(1, cols);

        // philipc: do need for dynamic list for Q
        //DynamicList<vectorField> Q(cols);
        List<vectorField> Q(cols);

        for (label i = 0; i < cols; i++)
        {
            //Q.append(QuasiNewtonV_[cols - 1 - i]);
            Q[i] = V_[cols - 1 - i];
        }

        for (label i = 0; i < cols; i++)
        {
            // Normalize column i
            R[i][i] = Foam::sqrt(gSum(Q[i] & Q[i]));
            Q[i] /= max(R[i][i], small);

            // Orthogonalize columns to the right of column i
            for (label j = i+1; j < cols; j++)
            {
                R[i][j] = gSum(Q[i] & Q[j]);
                Q[j] -= R[i][j]*Q[i];
            }

            // Project minus the residual vector on the Q
            C[i][0] =
                gSum
                (
                    Q[i]
                    & (
                        D.prevIter().primitiveField()
                      - D.primitiveField()
                    )
                );
        }

        // Solve the upper triangular system
        for (label j = 0; j < cols; j++)
        {
            Rcolsum[0][j] = 0.0;
            for (label i = 0; i < (j + 1); i++)
            {
                Rcolsum[0][j] += cmptMag(R[i][j]);
            }
        }
        scalar epsilon = 1.0E-10*max(Rcolsum);
        for (label i = 0; i < cols; i++)
        {
            if (cmptMag(R[i][i]) > epsilon)
            {
                for (label j = i + 1; j < cols; j++)
                {
                    R[i][j] /= R[i][i];
                }
                C[i][0] /= R[i][i];
                R[i][i] = 1.0;
            }
        }
        for (label j = (cols - 1); j >= 0; j--)
        {
            if (cmptMag(R[j][j]) > epsilon)
            {
                for (label i = 0; i < j; i++)
                {
                    C[i][0] -= C[j][0]*R[i][j];
                }
            }
            else
            {
                C[j][0] = 0.0;
            }
        }

        // Update D
        for (label i = 0; i < cols; i++)
        {
            D.primitiveFieldRef() += W_[i]*C[cols - 1 - i][0];
        }

        D.correctBoundaryConditions();
    }
    else
    {
        // Fixed under-relaxation during startup
        D.relax();
    }
}


// ************************************************************************* //
