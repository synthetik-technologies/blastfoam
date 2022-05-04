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
#include "fvc.H"
#include "distributedTriSurfaceMesh.H"
#include "volPointInterpolation.H"
#include "meshSizeObject.H"
#include "zeroGradientFvPatchFields.H"
#include "gaussGrad.H"

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

template<>
const char* Foam::NamedEnum<Foam::levelSetModel::truncation, 3>::names[] =
{
    "none",
    "cutOff",
    "tanh"
};

const Foam::NamedEnum<Foam::levelSetModel::levelSetFunc, 2>
    Foam::levelSetModel::levelSetFuncNames_;

const Foam::NamedEnum<Foam::levelSetModel::truncation, 3>
    Foam::levelSetModel::truncationNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::levelSetModel::levelSetModel
(
    volScalarField& alpha,
    const dictionary& dict,
    const bool mustRead
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
        dimensionedScalar("0", dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    H_
    (
        IOobject
        (
            IOobject::groupName("H", alpha.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    nHatf_
    (
        IOobject
        (
            IOobject::groupName("nHatf", alpha_.group()),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimArea, 0)
    ),
    K_
    (
        IOobject
        (
            IOobject::groupName("curvature", alpha_.group()),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    epsilon0_("epsilon", dimless, dict),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("levelSet::epsilon", alpha_.group()),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimLength, 0)
    ),
    useDistributed_(dict.lookupOrDefault("useDistributed", true)),
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
    ),
    truncation_
    (
        dict.found("truncation")
      ? truncationNames_.read(dict.lookup("truncation"))
      : truncation::NONE
    ),
    cutOff_
    (
        truncation_ == truncation::CUTOFF
      ? dict.lookup<scalar>("cutOffValue")
      : 0.0
    )
{
    updateEpsilon();
    if (Pstream::parRun())
    {
        triMeshDict_.set
        (
            "distributionType",
            distributedTriSurfaceMesh::distributionTypeNames_
            [
                distributedTriSurfaceMesh::FROZEN
            ]
        );
        triMeshDict_.set("mergeDistance", min(epsilon_).value()*1e-3);
    }

    if (levelSet_.headerOk())
    {
        correct(true);
    }
    else if (H_.headerOk())
    {
        correct(false);
    }
    else
    {
        if (mustRead)
        {
            levelSet_.readOpt() = IOobject::MUST_READ;
            levelSet_.read();
        }
        levelSet_ = calcLevelSet(alpha, 0.5);
        correct(true);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::levelSetModel::~levelSetModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::levelSetModel::updateEpsilon()
{
    epsilon_ = min(meshSizeObject::New(mesh_).dx())*epsilon0_;
}


Foam::tmp<Foam::volScalarField>
Foam::levelSetModel::calcH(const volScalarField& ls) const
{
    switch (lsFunc_)
    {
        case levelSetFunc::TANH:
        {
            return max(-1.0, min(1.0, tanh(ls/(2.0*epsilon_))));
        }
        case levelSetFunc::EXP:
        {
            return max(-1.0, min(1.0, 1.0 - 2.0/(1.0 + exp(ls/epsilon_))));
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown level set function" << endl
                << abort(FatalError);
        }
    }
    return ls;
}

Foam::tmp<Foam::volScalarField>
Foam::levelSetModel::calcLevelSet(const volScalarField& H) const
{
    switch (lsFunc_)
    {
        case levelSetFunc::TANH:
        {
            return atanh(max(small - 1.0, min(H, 1.0 - small)))*2.0*epsilon_;
        }
        case levelSetFunc::EXP:
        {
            return log(max(2.0/max(1.0 - H, small) - 1.0, small))*epsilon_;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown level set function" << endl
                << abort(FatalError);
        }
    }
    return H;
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::calcLevelSet
(
    const volScalarField& alpha,
    const UPtrList<searchableSurface>& regions
) const
{
    const fvMesh& mesh = alpha.mesh();
    tmp<volScalarField> tls
    (
        volScalarField::New
        (
            IOobject::groupName("levelSet", alpha.group()),
            mesh,
            dimensionedScalar("great", dimLength, great)
        )
    );
    volScalarField& ls = tls.ref();

    List<pointIndexHit> info(mesh.nCells());
    scalarField nearestDistSqr(mesh.nCells(), magSqr(mesh.bounds().span()));
    forAll(regions, regionI)
    {
        regions[regionI].findNearest(mesh.C(), nearestDistSqr, info);
        forAll(info, celli)
        {
            if (info[celli].hit())
            {
                ls[celli] =
                    min
                    (
                        ls[celli],
                        mag(info[celli].hitPoint() - mesh.C()[celli])
                    );
            }
        }
    }
    forAll(ls, celli)
    {
        if (alpha[celli] < 0.5)
        {
            ls[celli] *= -1.0;
        }
    }
    ls.correctBoundaryConditions();
    return tls;
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::calcLevelSet
(
    const volScalarField& isoField,
    const scalar isoValue
) const
{
    // Get point interpolated volume fraction field
    pointScalarField pointIsoField
    (
        volPointInterpolation::New(mesh_).interpolate(isoField)
    );

    // Contour the volume fraction at 0.5

    isoSurface contour
    (
        mesh_,
        isoField,
        pointIsoField,
        isoValue,
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
    if (Pstream::parRun() && useDistributed_)
    {
        triMeshDict_.set("bounds", List<boundBox>(1, mesh_.bounds()));
        triMeshPtr.set
        (
            new distributedTriSurfaceMesh
            (
                IOobject
                (
                    "contour_" + isoField.name(),
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
                    "contour_" + isoField.name(),
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
        triMesh.triSurface::write("contour_" + isoField.name() + ".stl");
    }

    // Temporary distance to interface field
    tmp<volScalarField> tls
    (
        volScalarField::New
        (
            IOobject::groupName("levelSet", isoField.group()),
            mesh_,
            dimensionedScalar(dimLength, -great),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& ls(tls.ref());

    // Collect all points that are need to be sampled (i.e. cell centers
    // and boundary face centres)
    pointField samples(mesh_.C());
    scalarField nearestDistSqr(samples.size(), magSqr(mesh_.bounds().span()));
    List<pointIndexHit> hitPoints(samples.size());

    //- Find the nearest points on the surface to the sample points
    triMesh.findNearest(samples, nearestDistSqr, hitPoints);

    // Compute distance to nearest point on the contour
    forAll(mesh_.C(), celli)
    {
        ls[celli] =
            mag(mesh_.C()[celli] - hitPoints[celli].rawPoint())
           *(alpha_[celli] > 0.5 ? 1.0 : -1.0);
    }
    ls.correctBoundaryConditions();

    return tls;
}


void Foam::levelSetModel::redistance()
{
    levelSet_ = calcLevelSet(levelSet_, 0.0);
}


void Foam::levelSetModel::correct(const bool updateH)
{
    updateEpsilon();

    if (updateH)
    {
        H_ = calcH(levelSet_);
    }
    else
    {
        levelSet_ = calcLevelSet(H_);
    }

    surfaceVectorField gradLevelSet(fvc::interpolate(fvc::grad(levelSet_)));
    nHatf_ =
        (
            gradLevelSet
           /max
            (
                mag(gradLevelSet),
                dimensionedScalar(gradLevelSet.dimensions(), 1e-6)
            )
        ) & mesh_.Sf();

    // Update curvature
    K_ = -fvc::div(nHatf_);
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::nearInterface() const
{
    return volScalarField::New
    (
        IOobject::groupName("nearInterface", levelSet_.group()),
        pos0(mag(levelSet_) + 2.0*epsilon_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::alpha() const
{
    tmp<volScalarField> tH(H_);
    switch (truncation_)
    {
        case truncation::NONE:
        {
            break;
        }
        case truncation::CUTOFF:
        {
            volScalarField cond(pos(mag(levelSet_) - sqrt(2.0)/2.0*epsilon_));
            tH = cond*sign(levelSet_) + (1.0 - cond)*levelSet_;
            break;
        }
        case truncation::TANH:
        {
            const scalar pi = Foam::constant::mathematical::pi;
            tH = tanh(pi*H_)/tanh(pi);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown truncation method" << endl
                << abort(FatalError);
        }
    }

    return volScalarField::New("alpha", 0.5*(tH + 1.0));
}

Foam::tmp<Foam::volVectorField> Foam::levelSetModel::gradAlpha() const
{
    // If truncation is used the gradient can be non-physical so force
    // linear
    if (truncation_ != truncation::NONE)
    {
        return fv::gaussGrad<scalar>(mesh_).grad(alpha());
    }
    return fvc::grad(alpha());
}


Foam::tmp<Foam::volVectorField> Foam::levelSetModel::nHat() const
{
    volVectorField gradLevelSet(fvc::grad(levelSet_));
    return volVectorField::New
    (
        "nHat",
        (
            gradLevelSet
           /max
            (
                mag(gradLevelSet),
                dimensionedScalar(gradLevelSet.dimensions(), 1e-6)
            )
        )
    );
}

// ************************************************************************* //
