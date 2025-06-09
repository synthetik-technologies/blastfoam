/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
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

#include "displacementConstraint.H"
#include "pointSet.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementConstraint, 0);
    defineRunTimeSelectionTable(displacementConstraint, dictionary);
}

template<>
const char* Foam::NamedEnum
<
    Foam::displacementConstraint::selectionType,
    3
>::names[] =
{
    "zone",
    "set",
    "patch"
};

const Foam::NamedEnum<Foam::displacementConstraint::selectionType, 3> Foam::displacementConstraint::selectionTypeNames;


Foam::List<Foam::label> Foam::displacementConstraint::readComponents
(
    Istream& is
)
{
    static HashTable<direction> validCmpts
    (
        {
            {"x", 0}, {"X", 0},
            {"y", 1}, {"Y", 1},
            {"z", 2}, {"Z", 2}
        }
    );
    List<word> cmptNames(is);
    labelHashSet cmpts;
    forAll(cmptNames, i)
    {
        HashTable<direction>::const_iterator iter =
            validCmpts.find(cmptNames[i]);
        if (iter == validCmpts.cend())
        {
            FatalIOErrorInFunction(is)
              <<  cmptNames[i] << " is not a valid component. " <<  nl
              << "Valid component names are: " << nl
              << validCmpts.sortedToc() << endl
              << abort(FatalIOError);
        }
        cmpts.insert(iter());
    }
    return cmpts.sortedToc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementConstraint::displacementConstraint
(
    const word& name,
    pointVectorField& D,
    pointVectorField& U,
    const dictionary& dict
)
:
    name_(name),
    mesh_(D.mesh().mesh()),
    D_(D),
    U_(U),
    selectionType_(UNKNOWN),
    selectionName_(name),
    dims_(readComponents(dict.lookup("dims"))),
    nodes_(0)
{
//     material id="1" name="dodgeball" type="neo-Hookean">
// 			<density>1000</density>
// 			<E>5516000</E>
// 			<v>0.32</v>
// 		</material>
// 		material id="2" name="steel cup" type="neo-Hookean">
// 			<density>7800</density>
// 			<E>2.1e+11</E>
// 			<v>0.3</v>
// 		</material>
//
// 		>0.001067
// 		10714
    // If selection was not specified and the name corresponds to a patch,
    // use the patch selection type
    if
    (
        !dict.found("selection")
     && mesh_.boundaryMesh().findIndex(name) >= 0
    )
    {
        selectionType_ = PATCH;
//         selectionName_ = name;
    }
    else
    {
        selectionType_ = selectionTypeNames.read(dict.lookup("selection"));
//         if (selectionType_ == ZONE)
//         {
//             selectionName_ = dict.lookup<word>("zone");
//         }
//         else if (selectionType_ == SET)
//         {
//             selectionName_ = dict.lookup<word>("set");
//         }
//         else if (selectionType_ == PATCH)
//         {
//             selectionName_ = dict.lookup<word>("patch");
//         }
//         else
//         {
//             FatalIOErrorInFunction(dict)
//                 << "Unknown selection type. Valid selectionTypes"
//                 << selectionTypeNames.toc() << endl
//                 << abort(FatalIOError);
//         }
    }
    updateNodes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementConstraint::~displacementConstraint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::displacementConstraint::updateNodes()
{
    switch (selectionType_)
    {
        case ZONE:
        {
            nodes_ = mesh_.pointZones()[selectionName_];
            return;
        }
        case SET:
        {
            nodes_ = pointSet(mesh_, selectionName_).toc();
            return;
        }
        case PATCH:
        {
            nodes_ = mesh_.boundaryMesh()
            [
                selectionName_
            ].meshPoints();
        }
        default:
        {
            // Need to handle
        }
    }
}
// ************************************************************************* //

