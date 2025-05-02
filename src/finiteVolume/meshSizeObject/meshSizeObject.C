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
    dXPtr_(nullptr),
    minDXPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshSizeObject::~meshSizeObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::meshSizeObject::movePoints()
{
    dxPtr_.clear();
    dXPtr_.clear();
    minDXPtr_.clear();
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

    dxPtr_.set(new scalarField(this->mesh().nCells(), 0.0));
    scalarField& dx = dxPtr_();
    const Vector<label>& geoD = this->mesh().geometricD();

    if (this->mesh().nGeometricD() == 1)
    {
        const labelListList& cellPoints = this->mesh().cellPoints();
        const pointField& points = this->mesh().points();
        label cmpti = -1;
        for (label i = 0; i < 3; i++)
        {
            if (geoD[i] > 0)
            {
                cmpti = i;
            }
        }
        forAll(cellPoints, celli)
        {
            dx[celli] = boundBox(points, cellPoints[celli], false).span()[cmpti];
        }
    }
    else if (this->mesh().nGeometricD() == 2)
    {
        forAll(this->mesh().boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchi];
            if (isA<wedgePolyPatch>(pp) || isA<emptyPolyPatch>(pp))
            {
                const List<label>& faceCells = pp.faceCells();
                forAll(faceCells, fi)
                {
                    dx[faceCells[fi]] += sqrt(pp.magFaceAreas()[fi]);
                }
            }
        }
        dx /= 2.0;
    }
    else
    {
        dx = cbrt(this->mesh().cellVolumes());
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
    dXPtr_.set(new vectorField(this->mesh().nCells(), vector::one));
    vectorField& dX = dXPtr_();

    const cellList& cells = this->mesh().cells();
    const scalarField& V = this->mesh().cellVolumes();
    const vectorField& Sf = this->mesh().faceAreas();

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


void Foam::meshSizeObject::calcMinDX() const
{
    if (minDXPtr_.valid())
    {
        FatalErrorInFunction
            <<"dX already set"
            << abort(FatalError);
    }
    minDXPtr_.set(new scalarField(this->mesh().nCells(), great));
    scalarField& minDX = minDXPtr_();
    const vectorField& DX = this->dX();

    for (label cmpti = 0; cmpti < 3; cmpti++)
    {
        if (this->mesh().geometricD()[cmpti] > 0)
        {
            forAll(minDX, celli)
            {
                minDX[celli] = min(minDX[celli], DX[celli][cmpti]);
            }
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::meshSizeObject::dx
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


Foam::tmp<Foam::volVectorField> Foam::meshSizeObject::dX
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
