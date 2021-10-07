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

#ifdef OPENFOAMESIORFOUNDATION

#include "extendedLeastSquaresVectors.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::extendedLeastSquaresVectors::extendedLeastSquaresVectors
(
    const fvMesh& mesh,
    const scalar minDet
)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, extendedLeastSquaresVectors>
    (
        mesh
    ),
    minDet_(minDet),
    pVectorsPtr_(nullptr),
    nVectorsPtr_(nullptr),
    additionalCellsPtr_(nullptr),
    additionalVectorsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::extendedLeastSquaresVectors::~extendedLeastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    deleteDemandDrivenData(additionalCellsPtr_);
    deleteDemandDrivenData(additionalVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extendedLeastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "extendedLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    pVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
#ifdef OPENFOAMESI
    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();
#else
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();
#endif
    const volVectorField& C = mesh().C();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh().nCells(), symmTensor::zero);

    forAll(owner, faceI)
    {
        // Build the d-vectors
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector d = C[nei] - C[own];

        symmTensor wdd = (1.0/magSqr(d))*sqr(d);

        dd[owner[faceI]] += wdd;
        dd[neighbour[faceI]] += wdd;
    }

    // Visit the boundaries. Coupled boundaries are taken into account
    // in the construction of d vectors.
    forAll(lsP.boundaryField(), patchI)
    {
        const fvPatch& p = mesh().boundary()[patchI];
#ifdef OPENFOAMESI
        const labelList& faceCells = p.faceCells();
#else
        const unallocLabelList& faceCells = p.faceCells();
#endif
        // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
        vectorField pd = p.delta();

        forAll(pd, patchFaceI)
        {
            dd[faceCells[patchFaceI]] +=
                (1.0/magSqr(pd[patchFaceI]))*sqr(pd[patchFaceI]);
        }
    }

    scalarField detdd = det(dd);

    label nAddCells = 0;
    label maxNaddCells = 4*detdd.size();
    additionalCellsPtr_ = new List<labelPair>(maxNaddCells);
    List<labelPair>& additionalCells_ = *additionalCellsPtr_;

    forAll (detdd, i)
    {
        label count = 0;

        while (++count < 100 && detdd[i] < minDet_)
        {
            if (nAddCells == maxNaddCells)
            {
                FatalErrorIn
                (
                    "extendedLeastSquaresVectors::"
                    "makeLeastSquaresVectors() const"
                )   << "nAddCells exceeds maxNaddCells"
                    << exit(FatalError);
            }

            labelList pointLabels = mesh().cells()[i].labels(mesh().faces());

            scalar maxDetddij = 0.0;

            label addCell = -1;

            forAll(pointLabels, j)
            {
                forAll(mesh().pointCells()[pointLabels[j]], k)
                {
                    label cellj = mesh().pointCells()[pointLabels[j]][k];

                    if (cellj != i)
                    {
                        if (cellj != -1)
                        {
                            vector dCij = (C[cellj] - C[i]);

                            symmTensor ddij =
                                dd[i] + (1.0/magSqr(dCij))*sqr(dCij);

                            scalar detddij = det(ddij);

                            if (detddij > maxDetddij)
                            {
                                addCell = cellj;
                                maxDetddij = detddij;
                            }
                        }
                    }
                }
            }

            if (addCell != -1)
            {
                additionalCells_[nAddCells][0] = i;
                additionalCells_[nAddCells++][1] = addCell;
                vector dCij = C[addCell] - C[i];
                dd[i] += (1.0/magSqr(dCij))*sqr(dCij);
                detdd[i] = det(dd[i]);
            }
        }
    }

    additionalCells_.setSize(nAddCells);

    // Invert the dd tensor
    symmTensorField invDd = inv(dd);


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, faceI)
    {
        // Build the d-vectors
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector d = C[nei] - C[own];

        lsP[faceI] = (1.0/magSqr(d))*(invDd[owner[faceI]] & d);

        lsN[faceI] =((-1.0)/magSqr(d))*(invDd[neighbour[faceI]] & d);
    }

    forAll(lsP.boundaryField(), patchI)
    {
        vectorField pd = mesh().boundary()[patchI].delta();

        fvsPatchVectorField& patchLsP = lsP.boundaryFieldRef()[patchI];

        const fvPatch& p = patchLsP.patch();
#ifdef OPENFOAMESI
        const labelList& faceCells = p.faceCells();
#else
        const unallocLabelList& faceCells = p.faceCells();
#endif

        forAll(p, patchFaceI)
        {
            patchLsP[patchFaceI] =
                (1.0/magSqr(pd[patchFaceI]))
               *(invDd[faceCells[patchFaceI]] & pd[patchFaceI]);
        }
    }


    additionalVectorsPtr_ = new vectorField(additionalCells_.size());
    vectorField& additionalVectors_ = *additionalVectorsPtr_;

    forAll(additionalCells_, i)
    {
        vector dCij = C[additionalCells_[i][1]]
            - C[additionalCells_[i][0]];

        additionalVectors_[i] =
            (1.0/magSqr(dCij))*(invDd[additionalCells_[i][0]] & dCij);
    }

    if (debug)
    {
        Info<< "extendedLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::surfaceVectorField&
Foam::extendedLeastSquaresVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField&
Foam::extendedLeastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


const Foam::List<Foam::labelPair>&
Foam::extendedLeastSquaresVectors::additionalCells() const
{
    if (!additionalCellsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *additionalCellsPtr_;
}


const Foam::vectorField&
Foam::extendedLeastSquaresVectors::additionalVectors() const
{
    if (!additionalVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *additionalVectorsPtr_;
}


bool Foam::extendedLeastSquaresVectors::movePoints()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    deleteDemandDrivenData(additionalCellsPtr_);
    deleteDemandDrivenData(additionalVectorsPtr_);

    return true;
}


#endif // end of #ifdef OPENFOAMESIORFOUNDATION

// ************************************************************************* //
