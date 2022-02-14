/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "levelSetModel.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "distributedTriSurfaceMesh.H"
#include "volPointInterpolation.H"
#include "meshSizeObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(levelSetModel, 0);
}

template<>
const char* Foam::NamedEnum<Foam::levelSetModel::levelSetFunc, 2>::names[] =
{
    "tanh",
    "exp"
};

const Foam::NamedEnum<Foam::levelSetModel::levelSetFunc, 2>
    Foam::levelSetModel::levelSetFuncNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::levelSetModel::levelSetModel
(
    const volScalarField& alpha,
    const dictionary& dict
)
:
    mesh_(alpha.mesh()),
    alpha_(alpha),
    levelSet_
    (
        IOobject
        (
            IOobject::groupName("levelSet", alpha.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        0.0
    ),
    epsilon_("epsilon", dimLength, dict),
    smoothPow_(dict.lookup<scalar>("smoothPow")),
    filterType_
    (
        dict.found("filtering")
      ? isoSurface::filterTypeNames_.read(dict.lookup("filtering"))
      : isoSurface::filterType::full
    ),
    lsFunc_
    (
        dict.found("levelSetFunction")
      ? levelSetFuncNames_.read(dict.lookup("levelSetFunction"))
      : levelSetFunc::TANH
    )
{
    scalar minDx(min(meshSizeObject::New(mesh_).dx()).value());
    if (epsilon_.value() < minDx)
    {
        WarningInFunction
            << "Interface thickness, " << epsilon_ << " is smaller than the "
            << "minimum mesh size, " << minDx << endl;
    }

    if (Pstream::parRun())
    {
        triMeshDict_.set
        (
            "distributionType",
            distributedTriSurfaceMesh::distributionTypeNames_
            [
                distributedTriSurfaceMesh::INDEPENDENT
            ]
        );
        triMeshDict_.set("mergeDistance", minDx/2.0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::levelSetModel::~levelSetModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::levelSetModel::calcLevelSet(const volScalarField& distance) const
{
    switch (lsFunc_)
    {
        case TANH:
        {
            return 0.5*(1.0 + tanh(distance/(2.0*epsilon_)));
        }
        case EXP:
        {
            return 1.0/(1.0 + exp(distance/epsilon_));
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown level set function" << endl
                << abort(FatalError);
        }
    }
}


void Foam::levelSetModel::redistance()
{
    // Get point interpolated volume fraction field
    pointScalarField pointAlpha
    (
        volPointInterpolation::New(mesh_).interpolate(alpha_)
    );

    // Contour the volume fraction at 0.5
    isoSurface contour
    (
        mesh_,
        alpha_,
        pointAlpha,
        0.5,
        filterType_
    );

    // Make sure the isoSurface is meshed with triangles
    contour.triangulate();

    // Copy faces to a triFaceList
    triFaceList triFaces(contour.size());
    forAll(contour, facei)
    {
        triFaces[facei][0] = contour[facei][0];
        triFaces[facei][1] = contour[facei][1];
        triFaces[facei][2] = contour[facei][2];
    }

    // Create a searchable triSufaceMesh
    autoPtr<triSurfaceMesh> triMeshPtr;
    triSurface tri(triFaces, contour.points());
    if (Pstream::parRun())
    {
        triMeshDict_.set("bounds", List<boundBox>(1, mesh_.bounds()));
        triMeshPtr.set
        (
            new distributedTriSurfaceMesh
            (
                IOobject
                (
                    IOobject::groupName("contour", alpha_.group()),
                    mesh_.time().timeName(),
                    mesh_
                ),
                tri,
                triMeshDict_
            )
        );
    }
    else
    {
        triMeshPtr.set
        (
            new triSurfaceMesh
            (
                IOobject
                (
                    IOobject::groupName("contour", alpha_.group()),
                    mesh_.time().timeName(),
                    mesh_
                ),
                tri
            )
        );
    }
    triSurfaceMesh& triMesh = triMeshPtr();

    if (debug)
    {
        triMesh.triSurface::write(alpha_.name() + "_contour.stl");
    }

    // Temporary distance to interface field
    volScalarField distance
    (
        volScalarField::New
        (
            IOobject::groupName("distance", alpha_.group()),
            mesh_,
            dimensionedScalar("0", dimLength, 0)
        )
    );

    // Calculate the distance to the cells from the nearest surface
    volScalarField::Boundary& bdistance = distance.boundaryFieldRef();

    // Collect all points that are need to be sampled (i.e. cell centers
    // and boundary face centres)
    pointField samples(mesh_.cellCentres());
    forAll(bdistance, patchi)
    {
        samples.append(mesh_.C().boundaryField()[patchi]);
    }
    scalarField nearestDistSqr(samples.size(), magSqr(mesh_.bounds().span()));
    List<pointIndexHit> hitPoints(samples.size());

    //- Find the nearest points on the surface to the sample points
    triMesh.findNearest(samples, nearestDistSqr, hitPoints);

    // Compute distance to nearest point on the contour
    forAll(mesh_.C(), celli)
    {
        distance[celli] =
            mag(mesh_.C()[celli] - hitPoints[celli].rawPoint())
           *(alpha_[celli] > 0.5 ? 1.0 : -1.0);
    }

    label i = mesh_.nCells();
    forAll(mesh_.boundaryMesh(), patchi)
    {
        scalarField& pdistance = bdistance[patchi];
        const vectorField& pCf = mesh_.Cf().boundaryField()[patchi];
        const scalarField& palpha = alpha_.boundaryField()[patchi];
        forAll(pdistance, facei)
        {
            pdistance[facei] =
                mag(pCf[facei] - hitPoints[i++].rawPoint())
               *(palpha[facei] > 0.5 ? 1.0 : -1.0);
        }
    }

    // Compute the level set function
    levelSet_ = calcLevelSet(distance);
    distance.write();
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::H() const
{
    return volScalarField::New
    (
        IOobject::groupName("H", alpha_.name()),
        tanh(sqr(alpha_*(1.0 - alpha_))/0.01)
    );
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::psi() const
{
    volScalarField lsPow(pow(levelSet_, smoothPow_));
    return volScalarField::New
    (
        IOobject::groupName("psi", alpha_.name()),
        lsPow/(lsPow + pow(1.0 - levelSet_, smoothPow_))
    );
}


Foam::tmp<Foam::volVectorField>
Foam::levelSetModel::gradLevelSet() const
{
    volScalarField alpha2(1.0 - alpha_);
    return volVectorField::New
    (
        "grad(" + IOobject::groupName("levelSet", alpha_.name() + ")"),
        fvc::grad(psi())/smoothPow_
       *pow(alpha_*alpha2, 1.0 - smoothPow_)
       *sqr(pow(alpha_, smoothPow_) + pow(alpha2, smoothPow_))
    );
}

// ************************************************************************* //
