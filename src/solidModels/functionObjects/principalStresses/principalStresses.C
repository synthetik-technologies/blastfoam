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
    defineTypeNameAndDebug(principalStresses, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        principalStresses,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::principalStresses::calculateEigenValues
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

bool Foam::principalStresses::writeData()
{
    if (runTime_.outputTime())
    {
        // Lookup stress tensor
        const volSymmTensorField& sigma =
            mesh_.lookupObject<volSymmTensorField>("sigma");

        // Calculate principal stress vectors

        volVectorField sigmaMax
        (
            IOobject
            (
                "sigmaMax",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("sigmaMaxVal", dimPressure, vector::zero)
        );

        volVectorField sigmaMid
        (
            IOobject
            (
                "sigmaMid",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("sigmaMaxVal", dimPressure, vector::zero)
        );

        volVectorField sigmaMin
        (
            IOobject
            (
                "sigmaMin",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("sigmaMaxVal", dimPressure, vector::zero)
        );

        // References to internalFields for efficiency
        const symmTensorField& sigmaI = sigma.primitiveField();
        vectorField& sigmaMaxI = sigmaMax.primitiveFieldRef();
        vectorField& sigmaMidI = sigmaMid.primitiveFieldRef();
        vectorField& sigmaMinI = sigmaMin.primitiveFieldRef();

        forAll (sigmaI, cellI)
        {
            calculateEigenValues
            (
                sigmaI[cellI],
                sigmaMaxI[cellI],
                sigmaMidI[cellI],
                sigmaMinI[cellI]
            );
        }

        forAll(sigmaMax.boundaryField(), patchI)
        {
            if
            (
                !sigmaMax.boundaryField()[patchI].coupled()
             && mesh_.boundaryMesh()[patchI].type() != "empty"
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
                }
            }
        }

        sigmaMax.correctBoundaryConditions();
        sigmaMid.correctBoundaryConditions();
        sigmaMin.correctBoundaryConditions();

        sigmaMax.write();
        sigmaMid.write();
        sigmaMin.write();

        Info<< "Principal stresses: max = " << gMax(mag(sigmaMax)()) << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::principalStresses::principalStresses
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    runTime_(t),
    mesh_
    (
        runTime_.lookupObject<fvMesh>
        (
            dict.lookupOrDefault<word>("region", "region0")
        )
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::principalStresses::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::principalStresses::execute()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::principalStresses::read(const dictionary& dict)
{
    return true;
}


bool Foam::principalStresses::write()
{
    return writeData();
}

// ************************************************************************* //
