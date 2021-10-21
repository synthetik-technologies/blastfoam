/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "meshSizeObject.H"
#include "fvc.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshSizeObject, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshSizeObject::meshSizeObject(const fvMesh& mesh)
:
    MeshSizeObject(mesh),
    mesh_(mesh),
    dxPtr_(nullptr),
    dXPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshSizeObject::~meshSizeObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshSizeObject::calcDx() const
{
    dxPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "dx",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimLength, 0.0)
        )
    );
    volScalarField& dx = dxPtr_();

    if (mesh_.nGeometricD() != 3)
    {

        scalarField faceLength(sqrt(mesh_.magFaceAreas()));
        const labelList& own = mesh_.faceOwner();
        labelList nFaces(dxPtr_->size(), 0);

        forAll(mesh_.boundaryMesh(), patchi)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchi];
            if (isA<wedgePolyPatch>(patch) || isA<emptyPolyPatch>(patch))
            {
                forAll(patch, fi)
                {
                    const label facei = patch.start() + fi;
                    dx[own[facei]] += faceLength[facei];

                    nFaces[own[facei]]++;
                }
            }
        }
        forAll(dx, celli)
        {
            dx[celli] /= scalar(nFaces[celli]);
        }
        forAll(mesh_.magSf().boundaryField(), patchi)
        {
            dx.boundaryFieldRef()[patchi] =
                dx.boundaryField()[patchi].patchInternalField();
        }
    }
    else
    {
        dx.primitiveFieldRef() = fvc::surfaceSum(mesh_.magSf())/mesh_.V();
    }

    forAll(mesh_.magSf().boundaryField(), patchi)
    {
        dx.boundaryFieldRef()[patchi] =
            sqrt(mesh_.magSf().boundaryField()[patchi]);
    }
}


void Foam::meshSizeObject::calcDX() const
{
    dXPtr_.set
    (
        new volVectorField
        (
            IOobject
            (
                "dX",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector(dimLength, Zero)
        )
    );
    volVectorField& dX = dXPtr_();

    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();
    const pointField& points = mesh_.points();
    forAll(dX, celli)
    {
        const cell& c = cells[celli];
        const edgeList cedges(c.edges(faces));
        vector dx(Zero);
        Vector<label> nEdges(Zero);
        forAll(cedges, ei)
        {
            const edge& e = cedges[ei];
            vector edx(cmptMag(e.vec(points)));
            label cmpti = findMax(edx);
            dx[cmpti] += edx[cmpti];
            nEdges[cmpti]++;
        }
        nEdges = max(nEdges, Vector<label>::one);
        dX[celli] = cmptDivide(dx, vector(nEdges));
    }

    for (label fi = mesh_.nInternalFaces(); fi < mesh_.nFaces(); fi++)
    {
        const label patchi = mesh_.boundaryMesh().whichPatch(fi);
        const polyPatch& p = mesh_.boundaryMesh()[patchi];
        const label facei = fi - p.start();

        const face& f = faces[fi];
        const edgeList& fedges = f.edges();
        vector dx(Zero);
        Vector<label> nEdges;
        forAll(fedges, ei)
        {
            const edge& e = fedges[ei];
            vector edx(cmptMag(e.vec(points)));
            label cmpti = findMax(edx);
            dx[cmpti] += edx[cmpti];
            nEdges[cmpti]++;
        }
        nEdges = max(nEdges, Vector<label>::one);
        dX.boundaryFieldRef()[patchi][facei] =
            cmptDivide(dx, vector(nEdges));
    }
}

// ************************************************************************* //
