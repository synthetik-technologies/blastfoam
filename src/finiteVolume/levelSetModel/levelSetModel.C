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
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha*2.0 - 1.0,
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
    smoothPow_(dict.lookupOrDefault<scalar>("smoothPow", 0.1)),
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
    if (truncation_ != truncation::NONE)
    {
        tlevelSet_.set
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("tlevelSet", alpha_.group()),
                    mesh_.time().timeName(),
                    mesh_
                ),
                levelSet_
            )
        );
    }

    scalar minDx(min(meshSizeObject::New(mesh_).dx()).value());
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
        triMeshDict_.set("mergeDistance", minDx/2.0);
    }
    correct(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::levelSetModel::~levelSetModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::levelSetModel::calcLevelSet(const volScalarField& d) const
{
    switch (lsFunc_)
    {
        case levelSetFunc::TANH:
        {
            return max(-1.0, min(1.0, tanh(d/(2.0*epsilon_))));
        }
        case levelSetFunc::EXP:
        {
            return max(-1.0, min(1.0, 1.0 - 2.0/(1.0 + exp(d/epsilon_))));
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown level set function" << endl
                << abort(FatalError);
        }
    }
}

Foam::tmp<Foam::volScalarField>
Foam::levelSetModel::calcDistance(const volScalarField& ls) const
{
    switch (lsFunc_)
    {
        case levelSetFunc::TANH:
        {
            return atanh(max(small - 1.0, min(ls, 1.0 - small)))*2.0*epsilon_;
        }
        case levelSetFunc::EXP:
        {
            return log(max(2.0/max(1.0 - ls, small) - 1.0, small))*epsilon_;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown level set function" << endl
                << abort(FatalError);
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::distance() const
{
    // Get point interpolated volume fraction field
    pointScalarField pointAlpha
    (
        volPointInterpolation::New(mesh_).interpolate(alpha_)
    );

    // Contour the volume fraction at 0.5
    contourPtr_.reset
    (
        new isoSurface
        (
            mesh_,
            alpha_,
            pointAlpha,
            0.5,
            filterType_
        )
    );
    isoSurface& contour = contourPtr_();

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
    tmp<volScalarField> td
    (
        volScalarField::New
        (
            IOobject::groupName("distance", alpha_.group()),
            mesh_,
            dimensionedScalar(dimLength, -great),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& d(td.ref());

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
        d[celli] =
            mag(mesh_.C()[celli] - hitPoints[celli].rawPoint())
           *(alpha_[celli] > 0.5 ? 1.0 : -1.0);
    }
    d.correctBoundaryConditions();

    return td;
}


const Foam::isoSurface& Foam::levelSetModel::contour() const
{
    if (!contourPtr_.valid())
    {
        distance();
    }
    return contourPtr_();
}


void Foam::levelSetModel::correct(const bool redistance)
{
    epsilon_ = min(meshSizeObject::New(mesh_).dx())*epsilon0_;
    if (redistance)
    {
        levelSet_ = calcLevelSet(distance());
    }
    if (tlevelSet_.valid())
    {
        switch (truncation_)
        {
            case truncation::NONE:
            {
                break;
            }
            case truncation::CUTOFF:
            {
                volScalarField d(calcDistance(levelSet_));
                volScalarField cond(pos(mag(d) - sqrt(2.0)/2.0*epsilon_));
                tlevelSet_() = cond*sign(levelSet_) + (1.0 - cond)*levelSet_;
                break;
            }
            case truncation::TANH:
            {
                const scalar pi = Foam::constant::mathematical::pi;
                tlevelSet_() = tanh(pi*levelSet_)/tanh(pi);
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown truncation method" << endl
                    << abort(FatalError);
            }
        }
    }


    surfaceVectorField gradPsi(fvc::interpolate(fvc::grad(psi())));
    nHatf_ =
        (
            gradPsi
           /max(mag(gradPsi), dimensionedScalar(gradPsi.dimensions(), 1e-6))
        ) & mesh_.Sf();

    // Update curvature
    K_ = -fvc::div(nHatf_);
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::nearInterface() const
{
    tmp<volScalarField> tnearInterface
    (
        volScalarField::New
        (
            IOobject::groupName("nearInterface", alpha_.name()),
            mesh_,
            0.0
        )
    );
    labelList intCells(interfaceCells());
    UIndirectList<scalar>(tnearInterface.ref(), intCells) = 1.0;
    return tnearInterface;
}


Foam::labelList Foam::levelSetModel::interfaceCells() const
{
    if (!contourPtr_.valid())
    {
        distance();
    }
    return contourPtr_->meshCells();
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::alpha() const
{
    return volScalarField::New("alpha", 0.5*(tlevelSet() + 1.0));
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::psi() const
{
    volScalarField alpha1(alpha());
    volScalarField alpha1Pow(pow(alpha1, smoothPow_));
    return volScalarField::New
    (
        IOobject::groupName("psi", alpha_.name()),
        alpha1Pow/(alpha1Pow + pow(1.0 - alpha1, smoothPow_))
    );
}


Foam::tmp<Foam::volVectorField>
Foam::levelSetModel::gradLevelSet() const
{
    volScalarField alpha1(0.5*(levelSet_ + 1.0));
    volScalarField alpha2(1.0 - alpha1);
    return volVectorField::New
    (
        "grad(" + IOobject::groupName("levelSet", alpha_.group() + ")"),
        fvc::grad(psi())/smoothPow_
       *pow(alpha1*alpha2, 1.0 - smoothPow_)
       *sqr(pow(alpha1, smoothPow_) + pow(alpha2, smoothPow_))
    );
}


Foam::tmp<Foam::volVectorField> Foam::levelSetModel::nHat() const
{
    volVectorField gradPsi(fvc::grad(psi()));
    return volVectorField::New
    (
        "nHat",
        (
            gradPsi
           /max(mag(gradPsi), dimensionedScalar(gradPsi.dimensions(), 1e-6))
        )
    );
}

// ************************************************************************* //
