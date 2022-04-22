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
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshSizeObject, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshSizeObject::meshSizeObject(const polyMesh& mesh)
:
    MeshSizeObject(mesh),
    dxPtr_(nullptr),
    dXPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshSizeObject::~meshSizeObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::meshSizeObject::movePoints()
{
    dxPtr_.clear();
    dXPtr_.clear();
    return true;
}


void Foam::meshSizeObject::calcDx() const
{
    if (dxPtr_.valid())
    {
        FatalErrorInFunction
            <<"dX already set"
            << abort(FatalError);
    }

    dxPtr_.set(new scalarField(mesh_.nCells()));
    scalarField& dx = dxPtr_();
    const Vector<label>& geoD = mesh_.geometricD();

    if (mesh_.nGeometricD() != 3)
    {
        const Vector<label>& solD = mesh_.solutionD();
        vector validD(Zero);
        forAll(solD, cmpti)
        {
            if (geoD[cmpti] < 0)
            {
                validD[cmpti] = 1.0;
            }
        }

        const vectorField& Sf = mesh_.faceAreas();
        const scalarField& magSf = mesh_.magFaceAreas();
        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();
        labelList nFaces(dxPtr_->size(), 0);

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            if (mag(Sf[facei]/magSf[facei] & validD) > 0.5)
            {
                dx[own[facei]] += magSf[facei];
                dx[nei[facei]] += magSf[facei];

                nFaces[own[facei]]++;
                nFaces[own[facei]]++;
            }
        }

        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            if (mag(Sf[facei]/magSf[facei] & validD) > 0.5)
            {
                dx[own[facei]] += magSf[facei];
                nFaces[own[facei]]++;
            }
        }

        forAll(dx, celli)
        {
            dx[celli] = sqrt(dx[celli]/scalar(nFaces[celli]));
        }
    }
    else
    {
        dx = cbrt(mesh_.cellVolumes());
    }
}


void Foam::meshSizeObject::calcDX() const
{
    if (dXPtr_.valid())
    {
        FatalErrorInFunction
            <<"dX already set"
            << abort(FatalError);
    }
    dXPtr_.set(new vectorField(mesh_.nCells(), vector::one));
    vectorField& dX = dXPtr_();

    const cellList& cells = mesh_.cells();
    const scalarField& V = mesh_.cellVolumes();
    const vectorField& Sf = mesh_.faceAreas();

    forAll(dX, celli)
    {
        const cell& c = cells[celli];
        vector sumMagSf(Zero);
        dX[celli] *= 2.0*V[celli];
        forAll(c, fi)
        {
            sumMagSf += cmptMag(Sf[c[fi]]);
        }
        dX[celli] = cmptDivide(dX[celli], sumMagSf);
    }
}


Foam::tmp<Foam::volScalarField> Foam::meshSizeObject::volDx
(
    const fvMesh& mesh
) const
{
    tmp<volScalarField> tdx
    (
        volScalarField::New
        (
            "dx",
            mesh,
            dimLength,
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    tdx.ref().primitiveFieldRef() = dx();
    tdx.ref().correctBoundaryConditions();
    return tdx;
}


Foam::tmp<Foam::volVectorField> Foam::meshSizeObject::volDX
(
    const fvMesh& mesh
) const
{
    tmp<volVectorField> tdX
    (
        volVectorField::New
        (
            "dX",
            mesh,
            dimLength,
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );
    tdX.ref().primitiveFieldRef() = dX();
    tdX.ref().correctBoundaryConditions();
    return tdX;
}


// ************************************************************************* //
