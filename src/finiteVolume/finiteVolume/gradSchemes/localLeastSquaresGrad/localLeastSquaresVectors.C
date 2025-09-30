/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "localLeastSquaresVectors.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(localLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::localLeastSquaresVectors::localLeastSquaresVectors
(
    const fvMesh& mesh,
    const bool local
)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        localLeastSquaresVectors
    >(mesh),
    local_(local),
    AinvPtr_(nullptr),
    AinvLocalPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::localLeastSquaresVectors::~localLeastSquaresVectors()
{
    deleteDemandDrivenData(AinvPtr_);
    deleteDemandDrivenData(AinvLocalPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::localLeastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        InfoInFunction
            << "Constructing least square gradient vectors"
            << endl;
    }

    AinvPtr_ = new tensorField(mesh().nCells(), Zero);
    tensorField& Ainv = *AinvPtr_;

    const volVectorField& C = mesh().C();
    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const vector dOwn = C[nei] - C[own];
        const vector dNei  = C[own] - C[nei];

        Ainv[own] += dOwn*dOwn;
        Ainv[nei] += dNei*dNei;
    }

    forAll(mesh().boundary(), patchi)
    {
        const fvPatch& patch = mesh().boundary()[patchi];
        const vectorField pd(patch.delta());

        if (patch.coupled())
        {
            forAll(patch, facei)
            {
                const label celli = patch.faceCells()[facei];
                Ainv[celli] += pd[facei]*pd[facei];
            }
        }
    }

    forAll(Ainv, celli)
    {
        Ainv[celli] = stabInv(Ainv[celli]);
    }

    if (debug)
    {
        InfoInFunction
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


void Foam::localLeastSquaresVectors::makeLocalLeastSquaresVectors() const
{
    if (debug)
    {
        InfoInFunction
            << "Constructing least square gradient vectors"
            << endl;
    }

    AinvPtr_ = new tensorField(mesh().nCells(), Zero);
    tensorField& Ainv = *AinvPtr_;

    const volVectorField& C = mesh().C();
    const surfaceVectorField& Cf = mesh().Cf();
    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const vector dOwn = Cf[facei] - C[own];
        const vector dNei = Cf[facei] - C[nei];

        Ainv[own] += dOwn*dOwn;
        Ainv[nei] += dNei*dNei;
    }

    forAll(mesh().boundary(), patchi)
    {
        const fvPatch& patch = mesh().boundary()[patchi];
        const vectorField pd(patch.fvPatch::delta());

        forAll(patch, facei)
        {
            const label celli = patch.faceCells()[facei];

            Ainv[celli] += pd[facei]*pd[facei];

            // if (D_.boundaryField()[patchi].fixesValue())
            // {
            //     const label mfacei =
            //         mesh_.boundary()[patchi].start() + facei;
            //     forAll(mesh_.faces()[mfacei], nodei)
            //     {
            //         const label pti = mesh_.faces()[mfacei][nodei];
            //
            //         vector d = points[pti] - C[celli];
            //         Ainv[celli] += d*d;
            //
            //         for (label i = 0; i < 7; i++)
            //         {
            //             scalar si(i);
            //             d =
            //                 (
            //                     (
            //                         (si + 1.0)*points[pti]
            //                       + (
            //                             (7.0 - si)
            //                            *Cf.boundaryField()[patchi][facei]
            //                         )
            //                     )
            //                 )/8.0 - C[celli];
            //             Ainv[celli] += d*d;
            //         }
            //     }
            // }
        }
    }
    forAll(Ainv, celli)
    {
        Ainv[celli] = stabInv(Ainv[celli]);
    }


    if (debug)
    {
        InfoInFunction
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::tensorField&
Foam::localLeastSquaresVectors::Ainv() const
{
    if (local_)
    {
        if (!AinvPtr_)
        {
            makeLeastSquaresVectors();
        }
        return *AinvPtr_;
    }

    if (!AinvPtr_)
    {
        makeLeastSquaresVectors();
    }
    return *AinvPtr_;

}


bool Foam::localLeastSquaresVectors::movePoints()
{
    deleteDemandDrivenData(AinvPtr_);
    deleteDemandDrivenData(AinvLocalPtr_);

    return true;
}


Foam::tensor Foam::localLeastSquaresVectors::stabInv(const tensor& t)
{
    if (magSqr(t) < small)
    {
        return tensor::zero;
    }

    scalar scale = magSqr(t);
    Vector<bool> removeCmpts
    (
        magSqr(t.xx())/scale < small,
        magSqr(t.yy())/scale < small,
        magSqr(t.zz())/scale < small
    );
    if (removeCmpts.x() || removeCmpts.y() || removeCmpts.z())
    {
        tensor tPlus(t);

        if (removeCmpts.x())
        {
            tPlus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tPlus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tPlus += tensor(0,0,0,0,0,0,0,0,1);
        }

        tensor tInv = inv(tPlus);

        if (removeCmpts.x())
        {
            tInv -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tInv -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tInv -= tensor(0,0,0,0,0,0,0,0,1);
        }
        return tInv;
    }
    return inv(t);
}

// ************************************************************************* //
