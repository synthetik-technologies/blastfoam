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

#include "hotCylinderAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hotCylinderAnalyticalSolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        hotCylinderAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::hotCylinderAnalyticalSolution::writeData()
{
    Info<< "hotCylinderAnalyticalSolution: writing analytical solution fields"
        << endl;

    // Lookup the solid mesh
    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    // Cell centre coordinates
    const volVectorField& C = mesh.C();

    // Create radial coordinates field
    // Note: I divide by 1.0 so that radii will be dimless
    volScalarField radii
    (
        sqrt
        (
            sqr(C.component(vector::X))
          + sqr(C.component(vector::Y))
        )/dimensionedScalar("one", dimLength, 1)
    );

    if (gMin(radii.primitiveField()) < SMALL)
    {
        FatalErrorIn("bool Foam::hotCylinderAnalyticalSolution::writeData()")
            << "The minimum pipe radius in zero: this is not allowed!"
            << " The pipe should be centered at the origin"
            << abort(FatalError);
    }

    // Create the analytical temperature field
    volScalarField analyticalT
    (
        IOobject
        (
            "analyticalT",
            time_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ((TInner_ - TOuter_)/Foam::log(rOuter_/rInner_))
       *Foam::log(rOuter_/radii)
    );

    // Write out the analytical temperature field
    Info<< "    Writing analytical temperature field (analyticalT)" << endl;
    analyticalT.write();

    // Create the analytical radial stress field
    volScalarField analyticalRadialStress
    (
        IOobject
        (
            "analyticalRadialStress",
            time_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (
            (alpha_*E_*(TInner_ - TOuter_))
           /(2.0*(1.0 - nu_)*Foam::log(rOuter_/rInner_))
        )
       *(
          - Foam::log(rOuter_/radii)
          - (
                sqr(rInner_)/(sqr(rOuter_) - sqr(rInner_))
            )*(1.0 - sqr(rOuter_)/sqr(radii))*Foam::log(rOuter_/rInner_)
        )
    );

    // Write out the analytical radial stress field
    Info<< "    Writing analytical radial stress field (analyticalRadialStress)"
        << endl;
    analyticalRadialStress.write();

    // Create the analytical hoop stress field
    volScalarField analyticalHoopStress
    (
        IOobject
        (
            "analyticalHoopStress",
            time_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (
            (alpha_*E_*(TInner_ - TOuter_))
           /(2.0*(1.0 - nu_)*Foam::log(rOuter_/rInner_))
        )
       *(
            1.0 - Foam::log(rOuter_/radii)
          - (
                sqr(rInner_)/(sqr(rOuter_) - sqr(rInner_))
            )*(1.0 + sqr(rOuter_)/sqr(radii))*Foam::log(rOuter_/rInner_)
        )
    );

    // Write out the analytical hoop stress field
    Info<< "    Writing analytical hoop stress field (analyticalHoopStress)"
        << nl << endl;
    analyticalHoopStress.write();

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hotCylinderAnalyticalSolution::hotCylinderAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    rInner_(readScalar(dict.lookup("rInner"))),
    rOuter_(readScalar(dict.lookup("rOuter"))),
    TInner_(readScalar(dict.lookup("TInner"))),
    TOuter_(readScalar(dict.lookup("TOuter"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu"))),
    alpha_(readScalar(dict.lookup("alpha")))
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (rInner_ >= rOuter_)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "rInner should be less than rOuter!"
            << abort(FatalError);
    }

    if (E_ < SMALL || nu_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "E and nu should be positive!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::hotCylinderAnalyticalSolution::start()
{
    return true;
}


bool Foam::hotCylinderAnalyticalSolution::execute()
{
    return writeData();
}


bool Foam::hotCylinderAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}


bool Foam::hotCylinderAnalyticalSolution::write()
{
    return writeData();
}

// ************************************************************************* //
