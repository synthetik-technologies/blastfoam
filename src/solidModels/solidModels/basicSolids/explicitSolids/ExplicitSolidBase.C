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

#include "ExplicitSolidBase.H"
#include "solidTractionFvPatchVectorField.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * * //

template<class IncrementalSolid>
void Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::relax()
{
    // Re-read the dictionary
    relaxation_.read
    (
        this->solidModelDict().optionalSubDict("relaxation")
    );

    //- Relax, if wanted
    relaxation_.relax(this->U(),this->rho());
    relaxation_.relax(a_);
    relaxation_.relax(this->DD());
    this->DD().correctBoundaryConditions();
}


template<class IncrementalSolid>
void Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::correctUBCs
(
    volVectorField& U
)
{
    U.correctBoundaryConditions();
    volVectorField::Boundary& bU = U.boundaryFieldRef();
    const volVectorField::Boundary& bDD = this->DD().boundaryField();
    forAll(bU, patchi)
    {
        if (isA<solidTractionFvPatchVectorField>(bDD[patchi]))
        {
            const fvPatch& patch = this->mesh().boundary()[patchi];

            fvPatchVectorField& pU = bU[patchi];
            const solidTractionFvPatchVectorField& pDD =
                dynamicCast<const solidTractionFvPatchVectorField>(bDD[patchi]);
            vectorField n(this->nf(patch));
            tensorField nn(n*n);

            tensorField St
            (
                nn/this->wavespeed_.boundaryField()[patchi]
            + (I - nn)/this->sWavespeed_.boundaryField()[patchi]
            );

            pU =
                pU.internalField()
              + (
                    St
                  & (
                        pDD.traction() - pDD.pressure()*n
                      - this->Pn(patch)
                    )
                )/this->rho().boundaryField()[patchi];
        }
        else if (bDD[patchi].fixesValue())
        {
            bU[patchi] = bDD[patchi]/this->mesh().time().deltaTValue();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalSolid>
Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::ExplicitSolidBase
(
    const word& type,
    fvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    IncrementalSolid(type, mesh, nonLinear, isSolid),
    wavespeed_
    (
        IOobject
        (
            "wavespeed",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVelocity, Zero)
    ),
    sWavespeed_
    (
        IOobject
        (
            "sWavespeed",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVelocity, Zero)
    ),
    energies_(mesh, this->solidModelDict()),
    a_
    (
        IOobject
        (
            "a",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity/dimTime, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    relaxation_(this->solidModelDict().optionalSubDict("relaxation"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalSolid>
void Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::eigenStructure
(
    const tensor& ten,
    vector& eigVal,
    tensor& eigVec
)
{
    tensor t(ten);

    scalar it_max = 100;
    scalar it_num = 0;
    scalar rot_num = 0;
    scalar size = 3;
    scalar gapq = 0.0;

    label i,j;
    label k,l,m,p1,q = 0;

    scalar term, termp, termq = 0.0;
    scalar g,h,c,w1 = 0.0;
    scalar theta,thresh = 0.0;
    scalar s,t1,tau1 = 0.0;

    vector d = vector(t.xx(),t.yy(),t.zz());
    vector bw = vector(t.xx(),t.yy(),t.zz());
    vector zw = vector::zero;
    tensor v1 = tensor::I;

    while (it_num++ < it_max)
    {
        // The convergence threshold is based on the size of the elements in
        // the strict upper triangle of the matrix.
        thresh = 0.0;

        for (j = 0; j < size; j++)
        {
            for (i = 0; i < j; i++)
            {
                thresh += sqr(t[i + j*size]);
            }
        }

        thresh = sqrt(thresh)/(4.0*size);

        if ( thresh < small )
        {
          break;
        }

        for (p1 = 0; p1 < size; p1++)
        {
            for (q = p1 + 1; q < size; q++)
            {
                gapq = 10.0*mag(t[p1 + q*size]);
                termp = gapq + mag(d[p1]);
                termq = gapq + mag(d[q]);

                // Annihilate tiny off-diagonal elements
                if (4 < it_num && termp == mag(d[p1]) && termq == mag(d[q]))
                {
                  t[p1+q*size] = 0.0;
                }

                //  Otherwise, apply a rotation
                else if (thresh <= mag(t[p1 + q*size]))
                {
                    h = d[q] - d[p1];
                    term = mag(h) + gapq;

                    if (term == mag(h))
                    {
                        t1 = t[p1 + q*size]/h;
                    }
                    else
                    {
                        theta = 0.5 * h/t[p1 + q*size];
                        t1 = 1.0/(mag(theta) + sqrt(1.0 + theta*theta));

                        if (theta < 0.0)
                        {
                            t1 = - t1;
                        }
                    }

                    c = 1.0/sqrt(1.0 + t1*t1);
                    s = t1*c;
                    tau1 = s/(1.0 + c);
                    h = t1*t[p1 + q*size];

                    //  Accumulate corrections to diagonal elements
                    zw[p1] = zw[p1] - h;
                    zw[q] = zw[q] + h;
                    d[p1] = d[p1] - h;
                    d[q] = d[q] + h;
                    t[p1 + q*size] = 0.0;

                    // Rotate, using information from the upper triangle of A
                    // only
                    for (j = 0; j < p1; j++)
                    {
                        g = t[j + p1*size];
                        h = t[j + q*size];
                        t[j + p1*size] = g - s*(h + g*tau1);
                        t[j + q*size] = h + s*(g - h*tau1);
                    }

                    for (j = p1 + 1; j < q; j++)
                    {
                        g = t[p1 + j*size];
                        h = t[j + q*size];
                        t[p1 + j*size] = g - s*(h + g*tau1);
                        t[j + q*size] = h + s*(g - h*tau1);
                    }

                    for (j = q + 1; j < size; j++)
                    {
                        g = t[p1 + j*size];
                        h = t[q + j*size];
                        t[p1 + j*size] = g - s*(h + g*tau1);
                        t[q + j*size] = h + s*(g - h*tau1);
                    }

                    //  Accumulate information in the eigenvector matrix
                    for (j = 0; j < size; j++)
                    {
                        g = v1[j + p1*size];
                        h = v1[j + q*size];
                        v1[j + p1*size] = g - s*(h + g*tau1);
                        v1[j + q*size] = h + s*(g - h*tau1);
                    }
                    rot_num = rot_num + 1;
                }
            }
        }
    }

    // Restore upper triangle of input matrix
    for (i = 0; i < size; i++)
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }

    //  Ascending sort the eigenvalues and eigenvectors
    for (k = 0; k < size - 1; k++)
    {
        m = k;
        for (l = k + 1; l < size; l++)
        {
            if (d[l] < d[m])
            {
                m = l;
            }
        }

        if (m != k)
        {
            t1 = d[m];
            d[m] = d[k];
            d[k] = t1;
            for (i = 0; i < size; i++)
            {
                w1 = v1[i + m*size];
                v1[i + m*size] = v1[i + k*size];
                v1[i + k*size] = w1;
            }
        }
    }

    // Corrections for calculating inverse
    tensor sub(tensor::zero);

    for (i = 0; i < 3; i++)
    {
        if (d[i] < SMALL)
        {
            d[i] = 1;
            sub +=
                d[i]
               *vector(v1[3*i], v1[3*i + 1], v1[3*i + 2])
               *vector(v1[3*i], v1[3*i + 1], v1[3*i + 2]);
        }
    }

    eigVal = d;
    eigVec = v1;
}


template<class IncrementalSolid>
void Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::solveMomentum()
{
    Info<< "Solving the momentum equation" << endl;

    this->enforceLinear(false);

    tmp<volVectorField> stab;


    // Central difference scheme
    const dimensionedScalar& deltaT = this->time().deltaT();
    const dimensionedScalar deltaT01
    (
        0.5*(deltaT + this->time().deltaT0())
    );

    // Compute the velocity
    // Note: this is the velocity at the middle of the time-step
    this->U() = this->U().oldTime() + deltaT01*a_.oldTime();

    // Compute change in displacement
    this->DD().internalFieldRef() = deltaT*this->U()();

    // Enforce any cell displacements
    if (this->setCellDisps().cellIDs().size())
    {
        vectorField& DDI = this->DD();
        vectorField& UI = this->U();
        const vectorField& DOld = this->D().oldTime();

        const labelList& cells = this->setCellDisps().cellIDs();
        const vectorField& cellDs = this->setCellDisps().cellDisps();

        forAll(cells, i)
        {
            const label celli = cells[i];
            DDI[celli] = cellDs[i] - DOld[celli];
            UI[celli] = DDI[celli]/deltaT.value();
        }
    }
    this->DD().correctBoundaryConditions();
    this->U().boundaryFieldRef() ==
        this->DD().boundaryField()/deltaT.value();
    // this->correctUBCs(this->U());

    relax();

    // Update displacement
    this->D() == this->D().oldTime() + this->U()*deltaT;

    // Update the stress field based on the latest D field
    this->update();

    // Compute acceleration
    // Note the inclusion of a linear bulk viscosity pressure term to
    // dissipate high frequency energies, and a Rhie-Chow term to
    // avoid checker-boarding
    stab =
    (
        this->stabilisation().stabilisation
        (
            this->U(),
            (deltaT01*this->impKf_)()
        )
    );
    a_ =
        (
            fvc::div(this->tractionSf())
          + fvc::div
            (
                this->mesh().Sf()*energies_.viscousPressuref
                (
                    this->rho(),
                    wavespeed_,
                    this->gradD()
                )
            )

            // This corresponds to Lax–Friedrichs smoothing
          + stab()
        )/this->rho()
      + this->g();
    a_.correctBoundaryConditions();


    // Check energies
    energies_.checkEnergies
    (
        this->rho(),
        this->U(),
        this->D(),
        this->DD(),
        this->sigma(),
        this->gradD(),
        this->gradDD(),
        this->stabilisation(),
        this->g()
    );
}


template<class IncrementalSolid>
Foam::scalar
Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::CoNum() const
{

    return
        this->mesh().time().deltaTValue()
       *gMax
        (
            this->mesh().surfaceInterpolation::deltaCoeffs().primitiveField()
           *this->wavespeed().primitiveField()
        );
/*
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value

    const fvMesh& mesh = this->mesh();
    surfaceScalarField amaxSf(this->wavespeed()*mesh.magSf());
    Info<<gMax(this->wavespeed())<<endl;

    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgeFvPatch>(mesh.boundary()[patchi]))
        {
            amaxSf.boundaryFieldRef() = Zero;
        }
    }
    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );
    return
        0.5*gMax(sumAmaxSf/mesh.V().field())*mesh.time().deltaTValue();*/
}


template<class IncrementalSolid>
Foam::scalar
Foam::solidModels::ExplicitSolidBase<IncrementalSolid>::maxCoNum() const
{
    return
        this->mesh().time().controlDict().lookupOrDefault
        (
            "maxCo",
            scalar(0.7071)
        );
}

// ************************************************************************* //
