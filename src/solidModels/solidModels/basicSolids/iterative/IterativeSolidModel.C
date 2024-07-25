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

template<class SolidModel>
void Foam::IterativeSolidModel<SolidModel>::relaxField
(
    volVectorField& D,
    const label iCorr
)
{
    if (relaxationMethod_ == "fixed")
    {
        // Fixed under-relaxation
        D.relax();
    }
    else if (relaxationMethod_ == "Aitken")
    {
        // See Aitken method at:
        // http://empire-multiphysics.com/projects/empire/wiki/Aitken_Relaxation
        // and
        // A partitioned solution approach for electro-thermo-
        // problems, Patrick Erbts, Stefan Hartmann, Alexander Duster.

        // Store aitkenResidual previous iteration
        aitkenResidual_().storePrevIter();

        // Calculate new aitkenResidual
        aitkenResidual_() = D.prevIter() - D;
        if (iCorr == 0)
        {
            // Fixed under-relaxation is applied on the first iteration
            aitkenAlpha_() = 1.0;

            if (this->mesh().relaxField(D.name()))
            {
                aitkenAlpha_() =
                    this->mesh().fieldRelaxationFactor(D.name());
            }
        }
        else
        {
            const volVectorField aitkenResidualDelta
            (
                aitkenResidual_().prevIter() - aitkenResidual_()
            );

            // Update the relaxation factor field
            aitkenAlpha_() =
                aitkenAlpha_()
               *(aitkenResidual_().prevIter() & aitkenResidualDelta)
               /(
                    magSqr(aitkenResidualDelta)
                  + dimensionedScalar("SMALL", dimLength*dimLength, SMALL)
                );

            // Bound alpha between 0.0 and 2.0
            // This may not be necessary but it seems to help convergence
            aitkenAlpha_() = max(0.0, min(2.0, aitkenAlpha_()));
        }

        // Relax the field
        D -= aitkenAlpha_()*aitkenResidual_();
    }
    else if (relaxationMethod_ == "QuasiNewton")
    {
        // This method is a modified form of the IQNILS by Degroote et al.

        // J. Degroote, K.-J. Bathe and J. Vierendeels.
        // A fluid solid interaction solver with IQN-ILS coupling algorithm.
        // Performance of a new partitioned procedure versus a monolithic
        // procedure in fluid-solid interaction. Computers & Solids

        if (iCorr == 0 || iCorr % QuasiNewtonRestartFreq_ == 0)
        {
            // Clean up data from old time steps

            if (SolidModel::debug)
            {
                Info<< "Modes before clean-up : " << QuasiNewtonT_.size();
            }

            while (true)
            {
                if (QuasiNewtonT_.size())
                {
                    if
                    (
                        this->runTime().timeIndex() > QuasiNewtonT_[0]
                     || iCorr % QuasiNewtonRestartFreq_ == 0
                    )
                    {
                        for (label i = 0; i < QuasiNewtonT_.size() - 1; i++)
                        {
                            QuasiNewtonT_[i] = QuasiNewtonT_[i + 1];
                            QuasiNewtonV_[i] = QuasiNewtonV_[i + 1];
                            QuasiNewtonW_[i] = QuasiNewtonW_[i + 1];
                        }

                        QuasiNewtonT_.remove();
                        QuasiNewtonV_.remove();
                        QuasiNewtonW_.remove();
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

            if (SolidModel::debug)
            {
                Info<< ", modes after clean-up : " << QuasiNewtonT_.size()
                    << endl;
            }
        }
        else if (iCorr == 1 || iCorr % QuasiNewtonRestartFreq_ == 1)
        {
            // Set reference in the first coupling iteration
            unrelaxedDRef_() = D;
            DRef_() = D.prevIter();
        }
        else
        {
            // Store the input vector field, defined as the previous iteration
            // D field (after relaxation) minus the Dp previous iteration field
            // in the first iteration (after relaxation)
            QuasiNewtonV_.append
            (
                (D - D.prevIter()) - (unrelaxedDRef_() - DRef_())
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
            QuasiNewtonW_.append(D - unrelaxedDRef_());

            // Store the time index
            QuasiNewtonT_.append(this->runTime().timeIndex());
        }

        if (QuasiNewtonT_.size() > 1)
        {
            // Consider QuasiNewtonV as a matrix V
            // with as columns the items
            // in the DynamicList and calculate the QR-decomposition of V
            // with modified Gram-Schmidt
            label cols = QuasiNewtonV_.size();
            RectangularMatrix<scalar> R(cols, cols, 0.0);
            RectangularMatrix<scalar> C(cols, 1);
            RectangularMatrix<scalar> Rcolsum(1, cols);
            // philipc: do need for dynamic list for Q
            //DynamicList<vectorField> Q(cols);
            List<vectorField> Q(cols);

            for (label i = 0; i < cols; i++)
            {
                //Q.append(QuasiNewtonV_[cols - 1 - i]);
                Q[i] = QuasiNewtonV_[cols - 1 - i];
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
                D.primitiveFieldRef() += QuasiNewtonW_[i]*C[cols - 1 - i][0];
            }

            D.correctBoundaryConditions();
        }
        else
        {
            // Fixed under-relaxation during startup
            D.relax();
        }
    }
    else
    {
        FatalErrorInFunction
            << "relaxationMethod '" << relaxationMethod_ << "' unknown!"
            << " Options are fixed, Aitken or QuasiNewton"
            << abort(FatalError);
    }
}


template<class SolidModel>
void Foam::IterativeSolidModel<SolidModel>::readDict()
{
    const dictionary& dict = this->solidModelDict();
    dict.readIfPresent("absTolerance", tolerance_);
    dict.readIfPresent("relTolerance", relTol_);
    dict.readIfPresent("solutionTolerance", solutionTol_);
    dict.readIfPresent("alternativeTolerance", alternativeTol_);
    dict.readIfPresent("materialTolerance", materialTol_);
    dict.readIfPresent("materialRelTolerance", materialRelTol_);
    dict.readIfPresent("infoFrequency", infoFrequency_);
    dict.readIfPresent("nCorrectors", nCorr_);
    dict.readIfPresent("minCorrectors", minCorr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidModel>
Foam::IterativeSolidModel<SolidModel>::IterativeSolidModel
(
    const word& type,
    dynamicFvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonlinear,
    const bool incremental,
    const bool isSolid
)
:
    SolidModel(type, mesh, nonlinear, incremental, isSolid),
    tolerance_(0.0),
    relTol_(1e-06),
    solutionTol_(1e-06),
    materialTol_(0.0),
    materialRelTol_(1e-05),
    alternativeTol_(1e-07),
    infoFrequency_(100),
    nCorr_(10000),
    minCorr_(1),
    maxIterReached_(0),
    residualFilePtr_(),
    relaxationMethod_
    (
        this->solidModelDict().template lookupOrDefault<word>
        (
            "relaxationMethod",
            "fixed"
        )
    ),
    aitkenAlpha_(),
    aitkenResidual_(),
    QuasiNewtonRestartFreq_(25),
    QuasiNewtonV_(),
    QuasiNewtonW_(),
    QuasiNewtonT_(),
    DRef_(),
    unrelaxedDRef_()
{
    readDict();

    // Print out the relaxation factor
    Info<< "    under-relaxation method: " << relaxationMethod_ << endl;
    if (relaxationMethod_ == "Aitken")
    {
        aitkenAlpha_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "aitkenAlpha",
                    mesh.time().constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("one", dimless, 1.0)
            )
        );
        aitkenResidual_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "aitkenResidual",
                    mesh.time().constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimLength, vector::zero)
            )
        );
    }
    else if (relaxationMethod_ == "QuasiNewton")
    {
        this->solidModelDict().readIfPresent
        (
            "QuasiNewtonRestartFrequency",
            QuasiNewtonRestartFreq_
        );
        Info<< "        restart frequency: " << QuasiNewtonRestartFreq_ << endl;
        QuasiNewtonV_.setSize(QuasiNewtonRestartFreq_ + 2),
        QuasiNewtonW_.setSize(QuasiNewtonRestartFreq_ + 2),
        QuasiNewtonT_.setSize(QuasiNewtonRestartFreq_ + 2),
        DRef_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "DRef",
                    mesh.time().constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimLength, vector::zero)
            )
        );
        unrelaxedDRef_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "unrelaxedDRef",
                    mesh.time().constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimLength, vector::zero)
            )
        );
    }

    // If requested, create the residual file
    if
    (
        this->solidModelDict().template lookupOrDefault<Switch>
        (
            "residualFile",
            false
        )
    )
    {
        if (Pstream::master())
        {
            Info<< "Creating residual.dat" << endl;
            residualFilePtr_.set
            (
                new OFstream(mesh.time().path()/"residual.dat")
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class SolidModel>
Foam::IterativeSolidModel<SolidModel>::~IterativeSolidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
