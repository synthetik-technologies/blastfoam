/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "surfacesToCell.H"
#include "boxToCell.H"
#include "cylinderToCell.H"
#include "sphereToCell.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfacesToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, surfacesToCell, word);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfacesToCell::surfacesToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    isBackup_(dict.dictName() == "backup")
{
    if (isBackup_)
    {
        if (dict.found("centre"))
        {
            backup_.set(new sphereToCell(mesh, dict));
        }
        else if (dict.found("point1") || dict.found("p1"))
        {
            backup_.set(new cylinderToCell(mesh, dict));
        }
        else if (dict.found("box") || dict.found("boxes"))
        {
            backup_.set(new sphereToCell(mesh, dict));
        }
        else
        {
            FatalErrorInFunction
                << "Unknown backup shape." << nl
                << "Must provide:" << nl
                << "    centre and radius" << nl
                << "    p1, p2, and radius" << nl
                << "    box or boxes" << endl
                << abort(FatalError);
        }
    }
    else
    {
        if (dict.found("surfaces"))
        {
            PtrList<entry> surfaces(dict.lookup("surfaces"));
            triSurfaces_.setSize(surfaces.size());
            triSurfaceSearches_.setSize(surfaces.size());
            surfaces_.setSize(surfaces.size());
            forAll(surfaces, i)
            {
                const dictionary& surfaceDict = surfaces[i].dict();
                const fileName name = surfaces[i].keyword();
                triSurfaces_.set(i, new triSurface(name));
                triSurfaceSearches_.set(i, new triSurfaceSearch(triSurfaces_[i]));

                transformer t(createTransformer(surfaceDict));
                if (t.transforms() || t.translates())
                {
                    pointField points(triSurfaces_[i].points());
                    t.transformPosition(points, points);
                    triSurfaces_[i].movePoints(points);
                }
                surfaces_.set
                (
                    i,
                    new surfaceToCell
                    (
                        mesh,
                        name,
                        triSurfaces_[i],
                        triSurfaceSearches_[i],
                        surfaceDict.lookup("outsidePoints"),
                        surfaceDict.lookupOrDefault("includeCut", true),
                        surfaceDict.lookupOrDefault("includeInside", true),
                        surfaceDict.lookupOrDefault("includeOutside", false),
                        surfaceDict.lookupOrDefault("useSurfaceOrientation", false),
                        surfaceDict.lookupOrDefault("nearDist", 0.0),
                        surfaceDict.lookupOrDefault("curvature", 0.0)
                    )
                );
            }
        }
        else
        {
            const fileName name = dict.lookup<fileName>("file").expand();
            triSurfaces_.append(new triSurface(name));
            triSurfaceSearches_.append(new triSurfaceSearch(triSurfaces_[0]));

            transformer t(createTransformer(dict));
            if (t.transforms() || t.translates())
            {
                pointField points(triSurfaces_[0].points());
                t.transformPosition(points, points);
                triSurfaces_[0].movePoints(points);
            }
            surfaces_.append
            (
                new surfaceToCell
                (
                    mesh,
                    name,
                    triSurfaces_[0],
                    triSurfaceSearches_[0],
                    dict.lookup("outsidePoints"),
                    dict.lookupOrDefault("includeCut", true),
                    dict.lookupOrDefault("includeInside", true),
                    dict.lookupOrDefault("includeOutside", false),
                    dict.lookupOrDefault("useSurfaceOrientation", false),
                    dict.lookupOrDefault("nearDist", 0.0),
                    dict.lookupOrDefault("curvature", 0.0)
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfacesToCell::~surfacesToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfacesToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (isBackup_)
    {
        backup_->applyToSet(action, set);
    }
    else
    {
        forAll(surfaces_, i)
        {
            surfaces_[i].applyToSet(action, set);
        }
    }
}


Foam::transformer Foam::surfacesToCell::createTransformer
(
    const dictionary& dict
) const
{
    transformer transforms;

    // Scale the
    if (dict.found("scale"))
    {
        const token t(dict.lookup("scale"));
        if (t.isScalar())
        {
            const scalar s(t.scalarToken());
            transforms =
                transformer::scaling(diagTensor(s, s, s))
              & transforms;
        }
        else
        {
            const vector v(dict.lookup("scale"));
            transforms =
                transformer::scaling(diagTensor(v.x(), v.y(), v.z()))
              & transforms;
        }
    }

    // Rotate
    if (dict.found("rotate"))
    {
        Pair<vector> n1n2(dict.lookup("rotate"));

        n1n2[0] /= mag(n1n2[0]);
        n1n2[1] /= mag(n1n2[1]);

        transforms =
            transformer::rotation(rotationTensor(n1n2[0], n1n2[1]))
            & transforms;
    }
    else if (dict.found("Rx"))
    {
        const scalar a(dict.lookup<scalar>("Rx"));
        transforms = transformer::rotation(Rx(degToRad(a))) & transforms;
    }
    else if (dict.found("Ry"))
    {
        const scalar a(dict.lookup<scalar>("Ry"));
        transforms = transformer::rotation(Ry(degToRad(a))) & transforms;
    }
    else if (dict.found("Rz"))
    {
        const scalar a(dict.lookup<scalar>("Rz"));
        transforms = transformer::rotation(Rz(degToRad(a))) & transforms;
    }
    else if (dict.found("Ra"))
    {
        Tuple2<vector, scalar> va(dict.lookup("Ra"));
        const vector v(va.first());
        const scalar a(va.second());
        transforms = transformer::rotation(Ra(v, degToRad(a))) & transforms;
    }

    // Translate
    if (dict.found("translate"))
    {
        const vector v(dict.lookup("translate"));
        transforms = transformer::translation(v) & transforms;
    }
    return transforms;
}

// ************************************************************************* //
