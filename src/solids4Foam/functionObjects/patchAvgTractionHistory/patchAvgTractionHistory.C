/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "patchAvgTractionHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "polyPatchID.H"
#include "surfaceFields.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchAvgTractionHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchAvgTractionHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::patchAvgTractionHistory::writeData()
{
    Info << "Writing average traction history" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volSymmTensorField& sigma =
        mesh.lookupObject<volSymmTensorField>("sigma");

    const symmTensorField& patchSigma =
        sigma.boundaryField()[patchIndex_];

//     const vectorField& patchS = mesh.Sf().boundaryField()[patchIndex_];

    const scalarField& patchMagS = mesh.magSf().boundaryField()[patchIndex_];

//     vector totForce =
//         gSum(patchS & patchSigma);

//     const volTensorField& gradD =
//         mesh.lookupObject<volTensorField>("grad(D)");

    vector totForce =
        gSum
        (
            mesh.Sf().boundaryField()[patchIndex_] & patchSigma
        );

//     vector totForce =
//         gSum
//         (
//             mesh.Sf().boundaryField()[patchIndex_]
//           & (patchSigma & (I+gradD.boundaryField()[patchIndex_]))
//         );

    vector avgTraction = totForce/(gSum(patchMagS) - SMALL);

    Info << "Average traction at patch " << patchName_
        << ": " << avgTraction << "(" << totForce
        << ")" << "\n" << endl;

    if (Pstream::master())
    {
        historyFilePtr_()
            << mesh.time().value()
                << tab << avgTraction.x()
                << tab << avgTraction.y()
                << tab << avgTraction.z()
                << tab << totForce.x()
                << tab << totForce.y()
                << tab << totForce.z() << endl;
    }

    if (time_.outputTime())
    {
        OFstream file
        (
            time_.timePath()/"tensile-sigmayy-sigmaxy-radius.dat"
        );

        file.precision(6);

        const vectorField Cf = mesh.boundary()[patchIndex_].Cf();

        file << "x" << tab << "sigmayy" << tab << "sigmaxy" << endl;
        forAll(Cf, faceI)
        {
            file << Cf[faceI].x() << tab
//                 << Cf[faceI].y() << tab
//                 << Cf[faceI].z() << tab
                << patchSigma[faceI].yy() << tab
                << patchSigma[faceI].xy() << endl;
        }
    }


    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchAvgTractionHistory::patchAvgTractionHistory
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_(),
    patchName_(dict.lookup("patchName")),
    patchIndex_(-1)
{
    Info << "Creating functio object " << name_ << "\n" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    polyPatchID patch(patchName_, mesh.boundaryMesh());

    if (!patch.active())
    {
        FatalErrorIn("patchAvgTractionHistory::patchAvgTractionHistory()")
            << "Patch name " << patchName_ << " not found."
                << abort(FatalError);
    }

    patchIndex_ = patch.index();

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                mesh.time().timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            fileName file("avgTraction_" + patchName_ + ".dat");

            // Open new file at start up
            historyFilePtr_.reset(new OFstream(historyDir/file));

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time"
                        << tab << "avgTractionX"
                        << tab << "avgTractionY"
                        << tab << "avgTractionZ"
                        << tab << "totForceX"
                        << tab << "totForceY"
                        << tab << "totForceZ" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchAvgTractionHistory::start()
{
    return writeData();
}


bool Foam::patchAvgTractionHistory::execute()
{
    return writeData();
}


bool Foam::patchAvgTractionHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::patchAvgTractionHistory::write()
{
    return writeData();
}

// ************************************************************************* //
