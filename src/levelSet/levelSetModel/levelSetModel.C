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
#include "fvm.H"
#include "distributedTriSurfaceMesh.H"
#include "volPointInterpolation.H"
#include "meshSizeObject.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "gaussGrad.H"
#include "isoSurface.H"
#include "upwind.H"
#include "fluxSchemeBase.H"

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
    const surfaceScalarField& phi,
    const dictionary& dict,
    const bool mustRead
)
:
    levelSetModel(alpha, dict, mustRead)
{
    phiPtr_.set(&phi);
}


Foam::levelSetModel::levelSetModel
(
    volScalarField& alpha,
    const dictionary& dict,
    const bool mustRead
)
:
    mesh_(alpha.mesh()),
    alpha_(alpha),
    phiPtr_(nullptr),
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
    epsilon_("epsilon", dimLength, 0.0),
    useDistributed_(dict.lookupOrDefault("useDistributed", false)),
    filterType_
    (
        dict.found("filtering")
      ? dict.lookup<word>("filtering")
      : isoSurface::filterTypeNames_[isoSurface::filterType::full]
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
    ),
    solveH_(true)
{
    volScalarField::Boundary& bls(levelSet_.boundaryFieldRef());
    forAll(bls, patchi)
    {
        if(isA<fixedGradientFvPatchScalarField>(bls[patchi]))
        {
            dynamicCast<fixedGradientFvPatchScalarField>
            (
                bls[patchi]
            ).gradient() = -1.0;
        }
    }

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
        triMeshDict_.set("mergeDistance", epsilon_.value()*1e-3);
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
    epsilon_ =
        dimensionedScalar
        (
            "dx",
            dimLength,
            min(meshSizeObject::New(mesh_).dx())
        )*epsilon0_;
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
        isoSurface::filterTypeNames_[filterType_]
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
           *(isoField[celli] > isoValue ? 1.0 : -1.0);
    }

    volScalarField::Boundary& bls = ls.boundaryFieldRef();
    forAll(bls, patchi)
    {
        const pointField& pCf = mesh_.Cf().boundaryField()[patchi];
        const scalarField& palpha = isoField.boundaryField()[patchi];
        triMesh.findNearest(pCf, nearestDistSqr, hitPoints);
        forAll(bls[patchi], facei)
        {
            bls[patchi][facei] =
                mag(pCf[facei] - hitPoints[facei].rawPoint())
               *(palpha[facei] > isoValue ? 1.0 : -1.0);
        }
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
    if (updateH)
    {
        H_ = calcH(levelSet_);
        H_.correctBoundaryConditions();
    }
    else
    {
        levelSet_ = calcLevelSet(H_);
        levelSet_.correctBoundaryConditions();
    }

    surfaceVectorField gradLevelSetf(fvc::interpolate(fvc::grad(levelSet_)));
    nHatf_ = (gradLevelSetf/max(mag(gradLevelSetf), small)) & mesh_.Sf();

    // Update curvature
    K_ = -fvc::div(nHatf_);
}


Foam::tmp<Foam::volScalarField> Foam::levelSetModel::nearInterface() const
{
    return volScalarField::New
    (
        IOobject::groupName("nearInterface", levelSet_.group()),
        pos0(epsilon_ - mag(levelSet_))
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
            volScalarField cond(pos(mag(H_) - cutOff_));
            tH = cond*sign(H_) + (1.0 - cond)*H_;
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


void Foam::levelSetModel::update()
{
    updateEpsilon();
}


void Foam::levelSetModel::solve()
{
    // const surfaceScalarField& phi = phiPtr_();
    // if (solveH_)
    // {
    //     tmp<surfaceScalarField> tHf;
    //     if (fluxSchemeBase::foundFluxScheme(phi))
    //     {
    //         tHf = fluxSchemeBase::findFluxScheme(phi).interpolate(H_);
    //     }
    //     else
    //     {
    //         tHf = fvc::interpolate(H_);
    //     }
    //     const surfaceScalarField& Hf(tHf());
    //
    //     surfaceScalarField gamma
    //     (
    //         IOobject::groupName("gamma", phi.group()),
    //         mag
    //         (
    //             phi/mesh_.magSf()
    //             // H_.mesh().lookupObject<volVectorField>
    //             // (
    //             //     IOobject::groupName("U", phi.group())
    //             // )
    //         )
    //     );
    //
    //
    //     volScalarField deltaH
    //     (
    //         IOobject::groupName("deltaH", alpha_.group()),
    //         fvc::div(phi*Hf)- H_*fvc::div(phi)
    //       // - fvc::laplacian(gamma*fvc::interpolate(epsilon_), H_)
    //       // + fvc::div(gamma*Hf*(1.0 - Hf)*nHatf_)
    //     );
    //
    //     volScalarField HOld(H_);
    //     this->storeAndBlendOld(HOld, false);
    //     this->storeAndBlendDelta(deltaH);
    //
    //     H_ = HOld - mesh_.time().deltaT()*deltaH;
    //     H_.maxMin(-1.0, 1.0);
    //     H_.correctBoundaryConditions();
    //
    //     deltaH = (H_ - HOld)/mesh_.time().deltaT();
    //     deltaH = this->calcAndStoreDelta(deltaH);
    //     correct(false);
    // }
    // else
    // {
    //     // Level set function is smooth so we do not need to worry about shocks
    //     // volScalarField deltaLevelSet
    //     // (
    //     //     IOobject::groupName("deltaLevelSet", alpha_.group()),
    //     //     fvc::div(phi, levelSet_)
    //     //   - levelSet_*fvc::div(phi)
    //     // );
    //     //
    //     // this->storeAndBlendOld(levelSet_, false);
    //     // this->storeAndBlendDelta(deltaLevelSet);
    //     // levelSet_ -= mesh_.time().deltaT()*deltaLevelSet;
    //     // levelSet_.correctBoundaryConditions();
    //     //
    //     // Info<<levelSet_.average().value()<<endl;
    //
    //     // H_ = calcH(levelSet_);
    //     // volScalarField levelSet(calcLevelSet(H_));
    //     dimensionedScalar dTau(dimLength, 0.1*min(meshSizeObject::New(mesh_).dx()));
    //     label nIter(epsilon_.value()/dTau.value());
    //     volScalarField S0("S0", sign(levelSet_));
    //     volScalarField nearInterface(this->nearInterface());
    //     if (mesh_.time().outputTime())
    //     {
    //         nearInterface.write();
    //     }
    //     for (label i = 0; i < nIter; i++)
    //     {
    //     //     levelSet.storePrevIter();
    //         levelSet_ += S0*(1.0 - mag(fvc::grad(levelSet_)))*dTau*nearInterface;
    //     }
    //     //     levelSet.correctBoundaryConditions();
    // //
    // //         // // levelSet_ += S0*(1.0 - mag(fvc::grad(levelSet_)));
    // //         //
    // //         // // volVectorField gradLevelSet(fvc::grad(levelSet_));
    // //         // // levelSet_ -= 0.5*S0*(1.0 - mag(gradLevelSet))/fvc::laplacian(S0, levelSet_);
    // //         //
    // //         // surfaceVectorField gradLevelSetf
    // //         // (
    // //         //     upwind<vector>(mesh_, fvc::interpolate(S0)).interpolate
    // //         //     (
    // //         //         fvc::grad(levelSet_)
    // //         //     )
    // //         // );
    // //         // volVectorField gradLevelSet(fvc::grad(levelSet_));
    // //         //
    // //         // surfaceScalarField w
    // //         // (
    // //         //     "w",
    // //         //     ((gradLevelSetf/max(mag(gradLevelSetf), small)) & mesh_.Sf())*S0f
    // //         // );
    // //         // volScalarField divW(fvc::div(w));
    // //         // // forAll(divW, celli)
    // //         // // {
    // //         // //     divW[celli] = stabilise(divW[celli], 1e-10);
    // //         // // }
    // //         //
    // //         // // levelSet_ -=
    // //         // //     (fvc::div(w, levelSet_) - fvc::div(w)*levelSet_ - S0)
    // //         // //    *dimensionedScalar(dimLength, 0.5);
    // //         // fvScalarMatrix levelSetEqn
    // //         // (
    // //         //     fvm::SuSp(volScalarField::New("one", mesh_, dx), levelSet_)
    // //         //   + fvm::div(w, levelSet_)
    // //         //   - fvm::Sp(divW, levelSet_)
    // //         //  ==
    // //         //     S0
    // //         // );
    // //         // levelSetEqn.relax(0.5);
    // //         // //
    // //         // //
    // //         // levelSetEqn.solve();
    // //         // // levelSet_.relax(0.5);
    // //         // Info<<levelSet_.average().value()<<" "<<(mag(levelSet_-levelSet_.prevIter()))().average().value()<<endl;
    // //         Info<<max(mag(mag(fvc::grad(levelSet))-1.0)).value()<<endl;
    // //     }
    //     // H_ = calcH(levelSet);
    // //
    //     correct(true);
    // }
    // // FatalErrorInFunctfion<<exit(FatalError);
    // // Info<<levelSet_.average()<<endl;
}


void Foam::levelSetModel::postUpdate()
{}


// ************************************************************************* //
