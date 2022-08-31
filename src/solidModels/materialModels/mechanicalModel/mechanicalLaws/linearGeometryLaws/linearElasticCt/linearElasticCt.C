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

#include "linearElasticCt.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "fvc.H"

#include "IFstream.H"
#include "bound.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticCt, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElasticCt, linGeomMechLaw
    );
}


// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

void Foam::linearElasticCt::setYoungsModulusFromCt()
{
    // Optional constant E
    if (dict().lookupOrDefault<Switch>("constantE", false))
    {
        Info<< nl << "Creating constant E field\n" << endl;

        E_ =
            dimensionedScalar
            (
                "E", dimPressure, readScalar(dict().lookup("value"))
            );

        Info<< nl << "Writing E field" << nl  << endl;

        E_.write();

        Info<< nl << "End" << endl;

    }
    else
    {
        Info<< nl << "Creating E field from CT Scans\n" << endl;

        // Lookup CT info from dictionary

        const int slices = readInt(dict().lookup("nSlices"));
        const int rows = readInt(dict().lookup("nRows"));
        const int columns = readInt(dict().lookup("nColumns"));
        const scalar pixelSliceSpacing =
            readScalar(dict().lookup("sliceSpacing"));
        const scalar pixelRowSpacing = readScalar(dict().lookup("rowSpacing"));
        const scalar pixelColumnSpacing =
            readScalar(dict().lookup("columnSpacing"));
        const scalar pixelSliceOffset =
            readScalar(dict().lookup("sliceOffset"));
        const scalar pixelRowOffset = readScalar(dict().lookup("rowOffset"));
        const scalar pixelColumnOffset =
            readScalar(dict().lookup("columnOffset"));

        //filepath location of CT image slices

        const fileName ctFilePath(dict().lookup("ctImagesFilePath"));

        Info<< "z: max " << slices*pixelSliceSpacing + pixelSliceOffset
            << " min " << pixelSliceOffset << nl
            << "y: max " << rows*pixelRowSpacing + pixelRowOffset
            << " min " << pixelRowOffset << nl
            << "x: max " << columns*pixelColumnSpacing + pixelColumnOffset
            << " min " << pixelColumnOffset << nl << endl;

        //- define ctImages array to hold ct HU values
        List<List<List<scalar> > > ctImages(slices); // use int
        forAll(ctImages, sliceI)
        {
            int percentComplete = 100*(sliceI+1)/slices;
            Info<< "\r\tDefining CT array to hold images: " << percentComplete
                << "%" << flush;

            ctImages[sliceI].setSize(rows);

            forAll(ctImages[sliceI], rowI)
            {
                ctImages[sliceI][rowI].setSize(columns, 0.0);
            }
        }

        Info<< nl << nl;

        Info<< "Reading in CT images" << endl;

        //- read ct files
        forAll(ctImages, sliceI)
        {
            int percentComplete = 100*(sliceI+1)/slices;
            Info<< "\r\tReading CT images: " << percentComplete << "%" << flush;
            std::stringstream ctFileName;

            ctFileName
                << ctFilePath <<"/ctImages/ct_scan_" << sliceI << ".txt";
            IFstream ctFile(ctFileName.str().c_str());

            if (!ctFile)
            {
                Info<< "Cannot open " << ctFileName.str().c_str() << "."
                    << endl;
            }

            forAll(ctImages[sliceI], rowI)
            {
                forAll(ctImages[sliceI][rowI], columnI)
                {
                    ctFile
                        >> ctImages[sliceI][rowI][columnI];
                }
            }
        }

        Info<< nl << "Creating Hu field\n" << endl;

        volScalarField Hu
        (
            IOobject
            (
                "Hu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("Hu0", dimPressure, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        Info<< nl << "Creating relRho field\n" << endl;

        volScalarField relRho
        (
            IOobject
            (
                "relRho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("relRho0", dimPressure, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        Info<< nl << "Creating E field\n" << endl;

        // Set the E field

        scalarField& E_I = E_.primitiveFieldRef();
        scalarField& HuI = Hu.primitiveFieldRef();
        scalarField& relRhoI = relRho.primitiveFieldRef();

        // Take a copy of the cell centres field as we may rotate them
        vectorField CI = mesh().C().primitiveField();

        if (useRotationMatrix_)
        {
            Info<< "    Rotating mesh coordinates before looking up CT values"
                << endl;
            CI =
                (rotationMatrix_ & (CI - centreOfRotation_))
              + centreOfRotation_;
        }

        forAll(E_I, cellI)
        {

            //-- Coordinates of the cell
            const scalar x_coord = CI[cellI].x();
            const scalar y_coord = CI[cellI].y();
            const scalar z_coord = CI[cellI].z();

            //C[cellI] = CI[cellI];

            //- set materialsI from CT
            int sliceNumber = 0;
            while
            (
                (z_coord >= (sliceNumber*pixelSliceSpacing + pixelSliceOffset))
             && sliceNumber < (slices-1)
            )
            {
                sliceNumber++;
            }

            // The slice number will now point to the slice just above the
            // z_coord
            // so the cell centre will be between (slice_number-1) and
            // (slice_number)

            //-- second find the two vertical pixel numbers the y_coord is
            // between
            int columnNumber=0;
            while
            (
                (
                     x_coord
                  >= (columnNumber*pixelColumnSpacing + pixelColumnOffset)
                )
             && columnNumber<(columns-1)
            )
            {
                columnNumber++;
            }
            // if I don't do this then it is the mirror image
            columnNumber = (columns-1) - columnNumber;

            // The vertical_pixel will now point to the pixal row just below
            // the of
            // y_coord so the y_coord will be between (vertical_pixel-1) and
            // (vertical_pixal)

            //-- second find the two horizontal pixel numbers the x_coord is
            // between
            int rowNumber = 0;
            while
            (
                (y_coord >= (rowNumber*pixelRowSpacing + pixelRowOffset))
             && rowNumber<(rows-1)
            )
            {
                rowNumber++;
            }

            // hmmnn I am not sure if I have to do this
            rowNumber = (rows - 1) - rowNumber;

            if
            (
               sliceNumber == 0 || rowNumber == 0 || columnNumber == 0
            || sliceNumber == (slices-1) || rowNumber == (rows-1)
            || columnNumber == (columns-1)
            )
            {
                HuI[cellI] = 0;
                relRhoI[cellI] = 0;
                E_I[cellI] = 0;
            }
            else
            {
                HuI[cellI] =
                    ctImages[sliceNumber][rowNumber][columnNumber];

                if ( HuI[cellI] > 816 )
                {
                    // Set function here
                    relRhoI[cellI] =
                        0.000769
                        *(ctImages[sliceNumber][rowNumber][columnNumber])
                      + 1.028;
                }
                else
                {
                    // Set function here
                    relRhoI[cellI] =
                        max
                        (
                            0.0019
                           *(ctImages[sliceNumber][rowNumber][columnNumber])
                          + 0.105,
                            0.2
                        );
                }

                if (relRhoI[cellI] > 1.54)
                {
                    // Set function here
                    E_I[cellI] = (0.09*(Foam::pow(relRhoI[cellI],7.4)))*1e9;
                }
                else
                {
                    // Set function here
                    E_I[cellI] =
                        (0.06 + 0.9*(Foam::pow(relRhoI[cellI], 2)))*1e9;
                }
            }
        }

        Hu.correctBoundaryConditions();
        E_.correctBoundaryConditions();
        relRho.correctBoundaryConditions();


        //Smooth the E field to help reduce oscillation due to the castellated
        //nature of the CT data

        // Optional smoothing
        if (dict().found("smooth"))
        {
            Info<< "Smoothing E field"<< endl;

            const scalar w = readScalar(dict().lookup("weight"));
            const label nSmoothIter = readInt(dict().lookup("nSmoothIters"));

            for (label i = 0; i < nSmoothIter; i++)
            {
                Info<< "    iteration " << i << endl;

                scalarField& E_I = E_.primitiveFieldRef();
                const vectorField& CI = mesh().C().primitiveField();
                const labelListList& cellCells = mesh().cellCells();

                forAll(E_I, cellI)
                {
                    const labelList& curCellCells = cellCells[cellI];

                    // Calculate average of neighbour cells
                    scalar av = 0.0;
                    scalar sumW = 0.0;
                    const vector& curCellC = CI[cellI];
                    forAll(curCellCells, ccI)
                    {
                        // Use inverse distance weighted average
                        scalar weight =
                            1.0/mag(curCellC - CI[curCellCells[ccI]]);
                        sumW += weight;

                        av += weight*E_I[curCellCells[ccI]];
                    }
                    av /= sumW;

                    // Set cell value to be weighted average of current cell
                    // value
                    // and that of the neighbours

                    E_I[cellI] = w*E_I[cellI] + (1.0 - w)*av;
                }
            }

            E_.correctBoundaryConditions();
        }

        // Optional: bound the E field
        if (dict().found("bound"))
        {
            Info<< "Bounding E field"<< endl;

            bound(E_, dimensionedScalar(dict().lookup("boundLowerValue")));
        }

        // Write the field

        Info<< nl << "Writing E field" << nl  << endl;

        E_.write();
    }

    Info<< nl << "End" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElasticCt::linearElasticCt
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    E_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar ("0.0", dimPressure, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    nu_("mu", dimless, dict),
    mu_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        E_/(1.0 + nu_)
    ),
    lambda_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        E_*nu_/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),
    muf_(fvc::interpolate(mu_)),
    lambdaf_(fvc::interpolate(lambda_)),
    useRotationMatrix_
    (
        dict.lookupOrDefault<Switch>("useRotationMatrix", false)
    ),
    rotationMatrix_(I),
    centreOfRotation_(vector::zero)
{
    if (useRotationMatrix_)
    {
        // Calculate rotation matrix using two vectors

        vector vectorBeforeRotation(dict.lookup("vectorBeforeRotation"));
        if (mag(vectorBeforeRotation) < SMALL)
        {
            FatalErrorIn(type() + "::" + type())
                << "mag(vectorBeforeRotation) < SMALL!" << abort(FatalError);
        }
        vectorBeforeRotation /= mag(vectorBeforeRotation);

        vector vectorAfterRotation(dict.lookup("vectorAfterRotation"));
        if (mag(vectorAfterRotation) < SMALL)
        {
            FatalErrorIn(type() + "::" + type())
                << "mag(vectorAfterRotation) < SMALL!" << abort(FatalError);
        }
        vectorAfterRotation /= mag(vectorAfterRotation);

        rotationMatrix_ =
            rotationTensor(vectorAfterRotation, vectorBeforeRotation);

        // Read centre of rotation
        centreOfRotation_ = vector(dict.lookup("centreOfRotation"));
    }

    // Set E field
    setYoungsModulusFromCt();

    // Reset mechanical fields after E has been updated
    mu_ = E_/(1.0 + nu_);
    lambda_ = E_*nu_/((1.0 + nu_)*(1.0 - 2.0*nu_));

    muf_ = fvc::interpolate(mu_);
    lambdaf_ = fvc::interpolate(lambda_);

    if (planeStress())
    {
        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));
        lambdaf_ = fvc::interpolate(lambda_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticCt::~linearElasticCt()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElasticCt::impK() const
{
    return volScalarField::New
    (
        "impK",
        2.0*mu_ + lambda_
    );
}


Foam::tmp<Foam::scalarField>
Foam::linearElasticCt::impK(const label patchi) const
{
    return
        2.0*mu_.boundaryField()[patchi]
      + lambda_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::linearElasticCt::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        lambda_ + 2.0/3.0*mu_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::linearElasticCt::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        lambda_ + 2.0*mu_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::linearElasticCt::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mu_
    );
}


void Foam::linearElasticCt::correct(volSymmTensorField& sigma)
{
    if (incremental())
    {
        FatalErrorIn
        (
            type() + "::correct(volSymmTensorField& sigma)"
        )   << "Not implemented for incremental solid solver"
            << abort(FatalError);
    }

    // Lookup gradient of displacment from the solver
    const volTensorField& gradD =
        mesh().lookupObject<volTensorField>("grad(D)");

    // Calculate stress
    sigma = mu_*twoSymm(gradD) + lambda_*tr(gradD)*I;
}


void Foam::linearElasticCt::correct(surfaceSymmTensorField& sigma)
{
    if (incremental())
    {
        FatalErrorIn
        (
            type() + "::correct(surfaceSymmTensorField& sigma)"
        )   << "Not implemented for incremental solid solver"
            << abort(FatalError);
    }

    // Lookup gradient of displacment from the solver
    const surfaceTensorField& gradD =
        mesh().lookupObject<surfaceTensorField>("grad(D)f");

    // Calculate stress
    sigma = muf_*twoSymm(gradD) + lambdaf_*tr(gradD)*I;
}



// ************************************************************************* //
