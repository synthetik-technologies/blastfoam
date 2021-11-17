/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "mechanicalModel.H"
#include "fvc.H"
#include "fvcGradf.H"
#include "gaussGrad.H"
#include "twoDPointCorrector.H"
#include "fixedGradientFvPatchFields.H"
#include "wedgePolyPatch.H"
#include "ZoneIDs.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mechanicalModel::makeVolToPoint() const
{
    if (volToPointPtr_)
    {
        FatalErrorIn
        (
            "void Foam::mechanicalModel::makeVolToPoint() const"
        )   << "pointer already set" << abort(FatalError);
    }

    volToPointPtr_ = new volPointInterpolation(mesh());
}


void Foam::mechanicalModel::calcImpKfcorr() const
{
    if (impKfcorrPtr_)
    {
        FatalErrorIn
        (
            "const Foam::volScalarField& "
            "Foam::mechanicalModel::calcImpKfcorr() const"
        )   << "pointer already set" << abort(FatalError);
    }

    impKfcorrPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "impKfcorr",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            impKf()
        );

    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() > 1)
    {
        FatalErrorIn(type())
            << "Not implemented for this version of OpenFOAM"
            << abort(FatalError);
    }
}


const Foam::surfaceScalarField& Foam::mechanicalModel::impKfcorr() const
{
    if (!impKfcorrPtr_)
    {
        calcImpKfcorr();
    }

    return *impKfcorrPtr_;
}


void Foam::mechanicalModel::clearOut()
{
    deleteDemandDrivenData(volToPointPtr_);
    deleteDemandDrivenData(impKfcorrPtr_);

    // Clear the list of mechanical laws
    // Note: we should do this before clearing the subMeshes, as the mechanical
    // laws can store geometricFields that must be deleted before deleting
    // mesh
    PtrList<mechanicalLaw>::clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalModel::mechanicalModel
(
    const fvMesh& mesh,
    const nonLinearGeometry::nonLinearType& nonLinGeom,
    const bool incremental
)
:
    IOdictionary
    (
        IOobject
        (
            "mechanicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    PtrList<mechanicalLaw>(),
    mesh_(mesh),
    planeStress_(lookup("planeStress")),
    incremental_(incremental),
    cellZoneNames_(),
    volToPointPtr_(),
    impKfcorrPtr_(NULL)
{
    Info<< "Creating the mechanicalModel" << endl;

    // Read the mechanical laws
    const PtrList<entry> lawEntries(lookup("mechanical"));

    PtrList<mechanicalLaw>& laws = *this;
    laws.setSize(lawEntries.size());

    // Create the list of cellZones names: they are used during the construction
    // of the subMeshes
    cellZoneNames_.setSize(laws.size());
    forAll(laws, lawI)
    {
        cellZoneNames_[lawI] = lawEntries[lawI].keyword();
    }

    // Create mechancial laws
    if (laws.size() == 1)
    {
        if (nonLinGeom == nonLinearGeometry::LINEAR_GEOMETRY)
        {
            laws.set
            (
                0,
                mechanicalLaw::NewLinGeomMechLaw
                (
                    lawEntries[0].keyword(),
                    mesh,
                    lawEntries[0].dict(),
                    nonLinGeom
                )
            );
        }
        else if
        (
            nonLinGeom == nonLinearGeometry::UPDATED_LAGRANGIAN
         || nonLinGeom == nonLinearGeometry::TOTAL_LAGRANGIAN
        )
        {
            laws.set
            (
                0,
                mechanicalLaw::NewNonLinGeomMechLaw
                (
                    lawEntries[0].keyword(),
                    mesh,
                    lawEntries[0].dict(),
                    nonLinGeom
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "It is not clear what type of mechanical law should be "
                << "created for a solidModel with nonLinGeom = " << nonLinGeom
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicalModel::~mechanicalModel()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::mechanicalModel::mesh() const
{
    return mesh_;
}


const Foam::volPointInterpolation& Foam::mechanicalModel::volToPoint() const
{
    if (!volToPointPtr_)
    {
        makeVolToPoint();
    }

    return *volToPointPtr_;
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::impK() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].impK();
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM"
            << abort(FatalError);

        return volScalarField::New
        (
            "impK",
            mesh_,
            dimForce/dimArea
        );
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalModel::impKf() const
{
    // Linear interpolation actually seems to give the best convergence
    const volScalarField impK(this->impK());
    const word interpName = "interpolate(" + impK.name() + ')';
    return fvc::interpolate(impK, interpName);
}


Foam::tmp<Foam::scalarField>
Foam::mechanicalModel::impK(const label patchi) const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].impK(patchi);
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM"
            << abort(FatalError);

        return tmp<scalarField>
        (
            new scalarField(mesh_.C().boundaryField()[patchi].size())
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::bulkModulus() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].bulkModulus();
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM"
            << abort(FatalError);

        return volScalarField::New
        (
            "bulkModulusLaw",
            mesh_,
            dimPressure
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::elasticModulus() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].elasticModulus();
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM"
            << abort(FatalError);

        return volScalarField::New
        (
            "elasticModulusLaw",
            mesh_,
            dimPressure
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalModel::shearModulus() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        return laws[0].shearModulus();
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM"
            << abort(FatalError);

        return volScalarField::New
        (
            "shearModulusLaw",
            mesh_,
            dimPressure
        );
    }
}


void Foam::mechanicalModel::correct(volSymmTensorField& sigma)
{
    PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        laws[0].correct(sigma);
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM" << abort(FatalError);
    }
}


void Foam::mechanicalModel::correct(surfaceSymmTensorField& sigma)
{
    PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        laws[0].correct(sigma);
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM" << abort(FatalError);
    }
}


void Foam::mechanicalModel::grad
(
    const volVectorField& D,
    volTensorField& gradD
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        gradD = fvc::grad(D);
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM" << abort(FatalError);
    }
}


void Foam::mechanicalModel::grad
(
    const volVectorField& D,
    const pointVectorField& pointD,
    volTensorField& gradD
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        gradD = fvc::grad(D, pointD);
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM" << abort(FatalError);
    }
}


void Foam::mechanicalModel::grad
(
    const volVectorField& D,
    const pointVectorField& pointD,
    surfaceTensorField& gradDf
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        gradDf = fvc::fGrad(D, pointD);
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM" << abort(FatalError);
    }
}


void Foam::mechanicalModel::grad
(
    const volVectorField& D,
    const pointVectorField& pointD,
    volTensorField& gradD,
    surfaceTensorField& gradDf
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        gradD = fvc::grad(D, pointD);
        gradDf = fvc::fGrad(D, pointD);
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM" << abort(FatalError);
    }
}


void Foam::mechanicalModel::interpolate
(
    const volVectorField& D,
    pointVectorField& pointD,
    const bool useVolFieldSigma
)
{
    const PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        volToPoint().interpolateDisplacement(D, pointD);
//         pointD.correctBoundaryConditions();
    }
    else
    {
        FatalErrorInFunction
            << "Not implemented for this version of OpenFOAM"
            << abort(FatalError);
    }
}


Foam::tmp<Foam::volVectorField> Foam::mechanicalModel::RhieChowCorrection
(
    const volVectorField& D,
    const volTensorField& gradD,
    const surfaceScalarField& gamma
) const
{
    // Mathematically "div(grad(phi))" is equivalent to "laplacian(phi)";
    // however, numerically "div(grad(phi))" uses a larger stencil than the
    // "laplacian(phi)"; the difference between these two approximations is
    // a small amount of numerical diffusion that quells oscillations
    //if (D.name() == "DD" || biMaterialInterfaceActive())
    if (true)
    {
        return
        (
            fvc::laplacian
            (
                gamma,
                D,
                "laplacian(D" + D.name() + ',' + D.name() + ')'
            )
          - fvc::div(gamma*mesh().Sf() & fvc::interpolate(gradD))
        );
    }
    else
    {
        // We will calculate this numerical diffusion based on the increment of
        // displacement, as it may become large of we base it on the total
        // displacement
        // Issue: The increment field "D - D.oldTime()" will be incorrect on
        // non-orthogonal meshes as the grad(D - D.oldTime()) field would not be
        // stored... we can/should fix this
        return
        (
            fvc::laplacian
            (
                gamma,
                D - D.oldTime(),
                "laplacian(D" + D.name() + ',' + D.name() + ')'
            )
          - fvc::div
            (
                gamma*mesh().Sf()
              & fvc::interpolate(gradD - gradD.oldTime())
            )
        );
    }
}


Foam::tmp<Foam::volVectorField> Foam::mechanicalModel::RhieChowCorrection
(
    const volVectorField& D,
    const volTensorField& gradD
) const
{
    return RhieChowCorrection(D, gradD, impKfcorr());
}


Foam::scalar Foam::mechanicalModel::residual()
{
    PtrList<mechanicalLaw>& laws = *this;

    scalar maxResidual = 0.0;

    forAll(laws, lawI)
    {
        maxResidual = max(maxResidual, laws[lawI].residual());
    }

    return maxResidual;
}


void Foam::mechanicalModel::updateTotalFields()
{
    PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        laws[lawI].updateTotalFields();
    }
}


Foam::scalar Foam::mechanicalModel::newDeltaT()
{
    // Find the minimum time-step of all the mechanical laws
    PtrList<mechanicalLaw>& laws = *this;

    // Initial set deltaT to as large as possible and then check
    // if any mechanical law wants a smaller time-step
    scalar newDeltaT = mesh().time().endTime().value();

    forAll(laws, lawI)
    {
        newDeltaT = min(newDeltaT, laws[lawI].newDeltaT());
    }

    return newDeltaT;
}


void Foam::mechanicalModel::moveSubMeshes()
{
}


void Foam::mechanicalModel::setUseSolidDeformation()
{
    PtrList<mechanicalLaw>& laws = *this;
    forAll(laws, lawI)
    {
        laws[lawI].setUseSolidDeformation();
    }
}

// ************************************************************************* //
