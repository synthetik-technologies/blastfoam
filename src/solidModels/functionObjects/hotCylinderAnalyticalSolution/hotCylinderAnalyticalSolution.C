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
namespace functionObjects
{
    defineTypeNameAndDebug(hotCylinderAnalyticalSolution, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        hotCylinderAnalyticalSolution,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::hotCylinderAnalyticalSolution::
hotCylinderAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, t, dict),
    rInner_("rInner", dimLength, dict.lookup("rInner")),
    rOuter_("rOuter", dimLength, dict.lookup("rOuter")),
    TInner_("TInner", dimTemperature, dict.lookup("TInner")),
    TOuter_("TOuter", dimTemperature, dict.lookup("TOuter")),
    E_("E", dimPressure, dict.lookup("E")),
    nu_("nu", dimless, dict.lookup("nu")),
    alpha_("alpha", dimVolume/dimTemperature, dict.lookup("alpha"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::hotCylinderAnalyticalSolution::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    dict.readIfPresent("rInner", rInner_.value());
    dict.readIfPresent("rOuter", rOuter_.value());
    dict.readIfPresent("TInner", TInner_.value());
    dict.readIfPresent("TOuter", TOuter_.value());
    dict.readIfPresent("E", E_.value());
    dict.readIfPresent("nu", nu_.value());
    dict.readIfPresent("alpha", alpha_.value());


    if (rInner_.value() >= rOuter_.value())
    {
        FatalErrorInFunction
            << "rInner should be less than rOuter!"
            << abort(FatalError);
    }

    if (E_.value() < small || nu_.value() < small)
    {
        FatalErrorInFunction
            << "E and nu should be positive!"
            << abort(FatalError);
    }

    return true;
}


bool Foam::functionObjects::hotCylinderAnalyticalSolution::execute()
{
    // Cell centre coordinates
    const volVectorField& C = mesh_.C();

    // Create radial coordinates field
    // Note: I divide by 1.0 so that radii will be dimless
    volScalarField radii
    (
        sqrt
        (
            sqr(C.component(vector::X))
          + sqr(C.component(vector::Y))
        )
    );

    if (gMin(radii.primitiveField()) < SMALL)
    {
        FatalErrorIn("bool Foam::hotCylinderAnalyticalSolution::writeData()")
            << "The minimum pipe radius in zero: this is not allowed!"
            << " The pipe should be centered at the origin"
            << abort(FatalError);
    }

    // Create the analytical temperature field
    store
    (
        "analyticalT",
        ((TInner_ - TOuter_)/Foam::log(rOuter_/rInner_))
       *Foam::log(rOuter_/radii)
    );

    // Create the analytical radial stress field
    store
    (
        "analyticalRadialStress",
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

    // Create the analytical hoop stress field
    store
    (
        "analyticalHoopStress",
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

    return true;
}


bool Foam::functionObjects::hotCylinderAnalyticalSolution::write()
{
    return
        writeObject("analyticalT")
     && writeObject("analyticalRadialStress")
     && writeObject("analyticalHoopStress");
}

// ************************************************************************* //
