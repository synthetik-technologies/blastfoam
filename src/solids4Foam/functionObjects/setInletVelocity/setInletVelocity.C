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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "setInletVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "boundBox.H"
#include "polyPatchID.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setInletVelocity, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setInletVelocity,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setInletVelocity::setVelocity()
{
    Info << "Seting inlet velocity" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    volVectorField& U =
        const_cast<volVectorField&>
        (
            mesh.lookupObject<volVectorField>("U")
        );


    word inletPatchName("inlet");

    polyPatchID inletPatch(inletPatchName, mesh.boundaryMesh());

    if (!inletPatch.active())
    {
        FatalErrorIn("setInletVelocity::setVelocity()")
            << "Inlet patch name " << inletPatchName << " not found."
            << abort(FatalError);
    }

    label inletPatchIndex = inletPatch.index();

    {
        fixedValueFvPatchVectorField& inletU =
            refCast<fixedValueFvPatchVectorField>
            (
                U.boundaryFieldRef()[inletPatchIndex]
            );

        const vectorField& Cf = mesh.Cf().boundaryField()[inletPatchIndex];

        scalar Umax = 0.2; //0.3;

        scalar t = time_.value() + time_.deltaT().value();

        scalar T = 1; //0.2;

        if (t < T)
        {
            Umax = Umax*(1.0 - cos(M_PI*t/T))/2.0;
        }

        Info << "Umax = " << Umax << endl;

        forAll(Cf, faceI)
        {
            scalar y = Cf[faceI].y();
            scalar z = Cf[faceI].z();

            inletU[faceI] =
                Umax*y*(0.4-y)*(sqr(0.4)-sqr(z))*vector(1,0,0)
               /(sqr(0.2)*sqr(0.4));
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setInletVelocity::setInletVelocity
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion)
{
    if (dict.found("region"))
    {
        dict.lookup("region")
            >> regionName_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setInletVelocity::start()
{
    return setVelocity();
}


bool Foam::setInletVelocity::execute()
{
    return setVelocity();
}


bool Foam::setInletVelocity::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region")
            >> regionName_;
    }

    return true;
}


bool Foam::setInletVelocity::write()
{
    return setVelocity();
}

// ************************************************************************* //
