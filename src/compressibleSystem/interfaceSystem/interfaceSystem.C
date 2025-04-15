/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "interfaceSystem.H"
#include "dynamicMultiAlphaContactAngleFvPatchScalarField.H"
#include "unitConversion.H"
#include "reconstruction.H"
#include "meshSizeObject.H"
#include "fvm.H"
#include "surfaceTensionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceSystem, 0);
    wordHashSet interfaceSystem::compressionSchemes
    (
        {
            "CICSAM",
            "HRIC",
            "mHRIC",
            "interfaceCompression",
            "PLIC",
            "PLICU",
            "MPLIC",
            "MPLICU"
        }
    );
}

// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::interfaceSystem::surfaceTensionForce
(
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const surfaceTensionModel& stModel
) const
{
    return
        surfaceScalarField::New
        (
            "surfaceTensionForce",
            fvc::interpolate(stModel.sigma()*K(alpha1, alpha2))*
            (
                fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
              - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
            )
        );
}


Foam::tmp<Foam::volVectorField> Foam::interfaceSystem::nHat
(
    const volScalarField& alpha
) const
{
    volVectorField gradAlpha(fvc::grad(alpha));

    // Unit interface normal
    return gradAlpha/(mag(gradAlpha) + deltaN_);
}


Foam::tmp<Foam::volVectorField> Foam::interfaceSystem::nHat
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    volVectorField gradAlpha
    (
        alpha2*fvc::grad(alpha1)
      - alpha1*fvc::grad(alpha2)
    );

    // Unit interface normal
    return gradAlpha/(mag(gradAlpha) + deltaN_);
}


Foam::tmp<Foam::surfaceVectorField> Foam::interfaceSystem::nHatfv
(
    const volScalarField& alpha
) const
{
    surfaceVectorField gradAlphaf(fvc::interpolate(fvc::grad(alpha)));

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceVectorField> Foam::interfaceSystem::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::interfaceSystem::nHatf
(
    const volScalarField& alpha
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha) & mesh_.Sf();
}


Foam::tmp<Foam::surfaceScalarField> Foam::interfaceSystem::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceSystem::correctContactAngle
(
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    surfaceVectorField::Boundary& nHatb
) const
{
    typedef dynamicMultiAlphaContactAngleFvPatchScalarField multiAlphaContact;
    const volScalarField::Boundary& a1bf = alpha1.boundaryField();
    const volScalarField::Boundary& a2bf = alpha2.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if
        (
            isA<multiAlphaContact>(a1bf[patchi])
         || isA<multiAlphaContact>(a2bf[patchi])
        )
        {
            if
            (
                isA<multiAlphaContact>(a1bf[patchi])
             && isA<multiAlphaContact>(a2bf[patchi])
            )
            {
                FatalErrorInFunction
                    << "alphaContactAngle boundary condition "
                       "specified on patch " << boundary[patchi].name()
                    << " for both " << alpha1.name() << " and " << alpha2.name()
                    << nl << "which may be inconsistent."
                    << exit(FatalError);
            }

            const multiAlphaContact& acap =
                isA<multiAlphaContact>(a1bf[patchi])
              ? refCast<const multiAlphaContact>(a1bf[patchi])
              : refCast<const multiAlphaContact>(a2bf[patchi])
              ;

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch
            (
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi]
            );

            multiAlphaContact::thetaPropsTable::
                const_iterator tp =
                acap.thetaProps().find(interfacePair(alpha1, alpha2));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            const bool matched = (tp.key().first() == alpha1.name());

            const scalar theta0 = degToRad(tp().theta0(matched));

            scalarField theta(boundary[patchi].size(), theta0);

            const scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > small)
            {
                const scalar thetaA = degToRad(tp().thetaA(matched));
                const scalar thetaR = degToRad(tp().thetaR(matched));

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    velocity_.boundaryField()[patchi].patchInternalField()
                  - velocity_.boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + small);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                const scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            const scalarField a12(nHatPatch & AfHatPatch);

            const scalarField b1(cos(theta));

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            const scalarField a((b1 - a12*b2)/det);
            const scalarField b((b2 - a12*b1)/det);

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
        else if (isA<alphaContactAngleFvPatchScalarField>(a1bf[patchi]))
        {
            FatalErrorInFunction
                << " When using more that 2 phases the "
                << multiAlphaContact::typeName
                << " boundary condition should be used in place of "
                << a1bf[patchi].type() << endl
                << abort(FatalError);
        }
        else if (isA<alphaContactAngleFvPatchScalarField>(a2bf[patchi]))
        {
            FatalErrorInFunction
                << " When using more that 2 phases the "
                << multiAlphaContact::typeName
                << " boundary condition should be used in place of "
                << a2bf[patchi].type() << endl
                << abort(FatalError);
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::interfaceSystem::K
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle(alpha1, alpha2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceSystem::nearInterface(const volScalarField& alpha) const
{
    return
        volScalarField::New
        (
            IOobject::groupName("nearInterface", alpha.group()),
            pos0(alpha - 0.01)*pos0(0.99 - alpha)
        );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceSystem::nearInterface(const surfaceScalarField& alphaf) const
{
    return
        surfaceScalarField::New
        (
            IOobject::groupName("nearInterface", alphaf.group()),
            pos0(alphaf - 0.01)*pos0(0.99 - alphaf)
        );
}


void Foam::interfaceSystem::updateDeltaN()
{
    deltaN_.value() = 1e-8/pow(average(mesh_.V()).value(), 1.0/3.0);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceSystem::nearInterface
(
    const volScalarField& psi,
    const dimensionedScalar& dX
) const
{
    dimensionedScalar dx(dimLength, average(meshSizeObject::New(mesh_).dx()));
    dimensionedScalar epsilon(1.5*dx);
    return
        volScalarField::New
        (
            IOobject::groupName("nearInterface", psi.group()),
            pos0(psi + epsilon)*pos0(psi - epsilon)
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceSystem::interfaceSystem
(
    const volVectorField& U,
    const dictionary& dict
)
:
    mesh_(U.mesh()),
    velocity_(U),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    ),
    dxs_({{"default", 0.005}})
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceSystem::~interfaceSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::interfaceSystem::usesCompression(const word& name) const
{
    word scheme
    (
        mesh_.interpolationScheme
        (
            reconstruction::scheme
            (
                IOobject::member(name),
                IOobject::group(name),
                mesh_,
                false,
                false
            )
        )
    );
    return compressionSchemes.found(scheme);
}


Foam::tmp<Foam::surfaceScalarField> Foam::interfaceSystem::surfaceTensionForce
(
    const UPtrList<volScalarField>& alphas,
    const surfaceTensionTable& stModels
) const
{
    tmp<surfaceScalarField> tstf
    (
        surfaceScalarField::New
        (
            "surfaceTensionForce",
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), 0.0)
        )
    );
    surfaceScalarField& stf = tstf.ref();

    forAll(alphas, phasei)
    {
        const volScalarField& alpha1 = alphas[phasei];
        forAll(alphas, phasej)
        {
            if (phasei == phasej)
            {
                continue;
            }
            const volScalarField& alpha2 = alphas[phasej];
            typename surfaceTensionTable::const_iterator iter =
                stModels.find(interfacePair(alpha1, alpha2));
            if (iter == stModels.cend())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of surface tension models"
                    << exit(FatalError);
            }

            stf +=
                fvc::interpolate(iter()->sigma()*K(alpha1, alpha2))
               *(
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }
    return tstf;
}


Foam::tmp<Foam::volScalarField> Foam::interfaceSystem::levelSet
(
    const volScalarField& alpha
) const
{
    dimensionedScalar dx(this->dx(alpha));
    dimensionedScalar gamma(0.75*dx);
    dimensionedScalar epsilon(1.5*dx);
    dimensionedScalar dTau(0.1*dx);

    tmp<volScalarField> tpsi
    (
        volScalarField::New
        (
            IOobject::groupName("levelSet", alpha.group()),
            (2.0*alpha - 1.0)*gamma,
            "zeroGradient"
        )
    );
    volScalarField& psi = tpsi.ref();

    const volScalarField S0(sign(psi));
    for (label i = 0; i < label(epsilon.value()/dTau.value()); i++)
    {
        psi += S0*(1.0 - mag(fvc::grad(psi)))*dTau;
        psi.correctBoundaryConditions();
    }

    return tpsi;
}


void Foam::interfaceSystem::redistance(volScalarField& psi) const
{
    dimensionedScalar dx(this->dx(psi));
    dimensionedScalar epsilon(1.5*dx);
    dimensionedScalar dTau(0.1*dx);

    const volScalarField S0(sign(psi));
    for (label i = 0; i < label(epsilon.value()/dTau.value()); i++)
    {
        psi += S0*(1.0 - mag(fvc::grad(psi)))*dTau;
        psi.correctBoundaryConditions();
    }
}


void Foam::interfaceSystem::GFMExtension
(
    volScalarField& x,
    const volScalarField& psi
) const
{
    // Gradient of level set function
    surfaceVectorField gradPsif(fvc::interpolate(fvc::grad(psi)));

    // Face unit interface normal
    surfaceScalarField nHatf((gradPsif & mesh_.Sf())/(mag(gradPsif) + 1e-6));
    volScalarField divnHatf(fvc::div(nHatf));

    GFMExtension(x, psi, nHatf, divnHatf);
}


void Foam::interfaceSystem::GFMExtension
(
    volScalarField& x,
    const volScalarField& psi,
    const surfaceScalarField& nHatf,
    const volScalarField& divnHatf
) const
{
    dimensionedScalar dx(this->dx(psi));
    dimensionedScalar dTau("dTau", 0.1*dx);

    const volScalarField xOld(x);
    for (label i = 0; i < 25; i++)
    {
        x -=
            dTau
            *(
                fvc::div(nHatf, x, "div(nHatf," + x.name() + ")")
              - x*divnHatf
            );

        x.correctBoundaryConditions();
    }

    x = x*pos0(psi) + xOld*neg(psi);
    x.correctBoundaryConditions();
}

// ************************************************************************* //
