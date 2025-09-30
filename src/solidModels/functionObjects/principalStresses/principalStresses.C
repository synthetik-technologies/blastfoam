/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "principalStresses.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(principalStresses, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        principalStresses,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::functionObjects::principalStresses::calculateEigenValues
(
    const symmTensor& sigma,
    vector& sigmaMax,
    vector& sigmaMid,
    vector& sigmaMin
)
{
    const vector eValues = eigenValues(sigma);
    const tensor eVectors = eigenVectors(sigma);

    label iMax = -1;
    label iMid = -1;
    label iMin = -1;
    label MaxMidMin = -1;
    // const label a = mag(eValues[0]);
    // const label b = mag(eValues[1]);
    // const label c = mag(eValues[2]);
    const scalar a = eValues[0];
    const scalar b = eValues[1];
    const scalar c = eValues[2];

    if (a < b)
    {
        if (a < c)
        {
            if (b < c)
            {
                iMin = 0;
                iMid = 1;
                iMax = 2;
            }
            else
            {
                iMin = 0;
                iMid = 2;
                iMax = 1;
            }
        }
        else
        {
            iMin = 2;
            iMid = 0;
            iMax = 1;
        }
    }
    else
    {
        if (b < c)
        {
            if (a < c)
            {
                iMin = 1;
                iMid = 0;
                iMax = 2;
            }
            else
            {
                iMin = 1;
                iMid = 2;
                iMax = 0;
            }
        }
        else
        {
            iMin = 2;
            iMid = 1;
            iMax = 0;
        }
    }

    MaxMidMin = iMax*100 + iMid*10 + iMin;

    if (MaxMidMin != -1)
    {
        switch (MaxMidMin)
        {
            case 12:
                sigmaMax = eVectors.x()*eValues.x();
                sigmaMid = eVectors.y()*eValues.y();
                sigmaMin = eVectors.z()*eValues.z();
                break;
            case 21:
                sigmaMax = eVectors.x()*eValues.x();
                sigmaMin = eVectors.y()*eValues.y();
                sigmaMid = eVectors.z()*eValues.z();
                break;
            case 102:
                sigmaMid = eVectors.x()*eValues.x();
                sigmaMax = eVectors.y()*eValues.y();
                sigmaMin = eVectors.z()*eValues.z();
                break;
            case 120:
                sigmaMid = eVectors.x()*eValues.x();
                sigmaMin = eVectors.y()*eValues.y();
                sigmaMax = eVectors.z()*eValues.z();
                break;
            case 201:
                sigmaMin = eVectors.x()*eValues.x();
                sigmaMax = eVectors.y()*eValues.y();
                sigmaMid = eVectors.z()*eValues.z();
                break;
            case 210:
                sigmaMin = eVectors.x()*eValues.x();
                sigmaMid = eVectors.y()*eValues.y();
                sigmaMax = eVectors.z()*eValues.z();
                break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::principalStresses::principalStresses
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, t, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::principalStresses::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}


bool Foam::functionObjects::principalStresses::execute()
{
    // Calculate principal stress vectors
    tmp<volVectorField> tsigmaMax
    (
        new volVectorField
        (
            IOobject
            (
                "sigmaMax",
                time_.name(),
                mesh_
            ),
            mesh_,
            dimensionedVector(dimPressure, vector::zero)
        )
    );
    volVectorField& sigmaMax = tsigmaMax.ref();


    tmp<volVectorField> tsigmaMin
    (
        new volVectorField
        (
            IOobject
            (
                "sigmaMin",
                time_.name(),
                mesh_
            ),
            mesh_,
            dimensionedVector(dimPressure, vector::zero)
        )
    );
    volVectorField& sigmaMin = tsigmaMin.ref();


    tmp<volVectorField> tsigmaMid
    (
        new volVectorField
        (
            IOobject
            (
                "sigmaMid",
                time_.name(),
                mesh_
            ),
            mesh_,
            dimensionedVector(dimPressure, vector::zero)
        )
    );
    volVectorField& sigmaMid = tsigmaMid.ref();

    // Lookup stress tensor
    const volSymmTensorField& sigma =
        mesh_.lookupObject<volSymmTensorField>("sigma");

    // References to internalFields for efficiency
    const symmTensorField& sigmaI = sigma.primitiveField();
    vectorField& sigmaMaxI = sigmaMax.primitiveFieldRef();
    vectorField& sigmaMidI = sigmaMid.primitiveFieldRef();
    vectorField& sigmaMinI = sigmaMin.primitiveFieldRef();

    scalar maxSigmaMag = 0.0;
    forAll (sigmaI, cellI)
    {
        calculateEigenValues
        (
            sigmaI[cellI],
            sigmaMaxI[cellI],
            sigmaMidI[cellI],
            sigmaMinI[cellI]
        );
        maxSigmaMag = max(maxSigmaMag, mag(sigmaMaxI[cellI]));
    }

    forAll(sigmaMax.boundaryField(), patchI)
    {
        if
        (
            !sigma.boundaryField()[patchI].coupled()
         && !isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchI])
        )
        {
            const symmTensorField& pSigma = sigma.boundaryField()[patchI];
            vectorField& pSigmaMax = sigmaMax.boundaryFieldRef()[patchI];
            vectorField& pSigmaMid = sigmaMid.boundaryFieldRef()[patchI];
            vectorField& pSigmaMin = sigmaMin.boundaryFieldRef()[patchI];

            forAll(pSigmaMax, faceI)
            {
                calculateEigenValues
                (
                    pSigma[faceI],
                    pSigmaMax[faceI],
                    pSigmaMid[faceI],
                    pSigmaMin[faceI]
                );

                maxSigmaMag = max(maxSigmaMag, mag(pSigmaMax[faceI]));
            }
        }
    }

    sigmaMax.correctBoundaryConditions();
    sigmaMid.correctBoundaryConditions();
    sigmaMin.correctBoundaryConditions();

    Info<< "Principal stresses: max = "
        << returnReduce(maxSigmaMag, maxOp<scalar>()) << endl;

    store(tsigmaMax);
    store(tsigmaMin);
    store(tsigmaMid);

    return true;
}


bool Foam::functionObjects::principalStresses::write()
{
    return
        writeObject("sigamMax")
     && writeObject("sigamMid")
     && writeObject("sigamMin");
}

// ************************************************************************* //
