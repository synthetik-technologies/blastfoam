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

#include "pointHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "IOmanip.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        pointHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::pointHistory::writeData()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    vector deflection = vector::zero;
    vector velocity = vector::zero;

    if (processor_ == Pstream::myProcNo())
    {
        deflection = mesh.points()[historyPointID_] - refHistoryPoint_;

        if (mesh.foundObject<pointVectorField>("pointD"))
        {
            const pointVectorField& pointU =
                mesh.lookupObject<pointVectorField>("pointD");

            deflection = pointU[historyPointID_];
        }
        else if (mesh.foundObject<pointVectorField>("pointU"))
        {
            const pointVectorField& pointU =
                mesh.lookupObject<pointVectorField>("pointU");

            deflection = pointU[historyPointID_];
        }
    }


    reduce(deflection, sumOp<vector>());
    reduce(velocity, sumOp<vector>());

    if (Pstream::master())
    {
        historyFilePtr_() << setprecision(12);

        historyFilePtr_()
            << time_.time().value() << tab
            << deflection.x() << tab
            << deflection.y() << tab
            << deflection.z();

        if (writeVelocity_)
        {
            historyFilePtr_()
                << tab << velocity.x() << tab
                    << velocity.y() << tab << velocity.z();
        }

        historyFilePtr_() << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointHistory::pointHistory
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
    historyPointID_(-1),
    refHistoryPoint_(dict.lookup("refHistoryPoint")),
    writeVelocity_(false),
    processor_(-1),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object." << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    if (dict.found("historyPointID"))
    {
        dict.lookup("historyPointID") >> historyPointID_;
    }

    if (dict.found("writeVelocity"))
    {
        dict.lookup("writeVelocity") >> writeVelocity_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const vectorField& points = mesh.points();

    List<scalar> minDist(Pstream::nProcs(), GREAT);

    if (historyPointID_ == -1)
    {
        forAll(points, pointI)
        {
            scalar dist = mag(refHistoryPoint_ - points[pointI]);

            if (dist < minDist[Pstream::myProcNo()])
            {
                minDist[Pstream::myProcNo()] = dist;
                historyPointID_ = pointI;
            }
        }
    }

    Pstream::gatherList(minDist);
    Pstream::scatterList(minDist);

    processor_ = -1;
    scalar min = GREAT;

    forAll(minDist, procI)
    {
        if (minDist[procI] < min)
        {
            min = minDist[procI];
            processor_ = procI;
        }
    }

    if (processor_ == Pstream::myProcNo())
    {
        Pout << "History point ID: " << historyPointID_ << endl;
        Pout << "History point coordinates: "
            << points[historyPointID_] << endl;
        Pout << "Reference point coordinates: " << refHistoryPoint_ << endl;
    }

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());


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



            // Open new file at start up

            OStringStream FileName;
            FileName() << "point_" << historyPointID_ << ".dat";

            historyFilePtr_.reset
            (
                new OFstream(historyDir/word(FileName.str()))
            );
//             historyFilePtr_.reset(new OFstream(historyDir/"point.dat"));

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << tab << "X" << tab
                        << "Y" << tab << "Z";

                if (writeVelocity_)
                {
                    historyFilePtr_()
                        << tab << "Vx" << tab
                            << "Vy" << tab << "Vz";
                }

                historyFilePtr_() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointHistory::start()
{
    return writeData();
}


bool Foam::pointHistory::execute()
{
    return writeData();
}


bool Foam::pointHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::pointHistory::write()
{
    return writeData();
}

// ************************************************************************* //
