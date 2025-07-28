/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
09-09-2024 Synthetik Applied Technologies: |    Calculate Courant number with
                                                blastFoam thermo
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

#include "blastCourantNo.H"
#include "fluidBlastThermo.H"
#include "fvcMeshPhi.H"
#include "fvc.H"
#include "wedgePolyPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(blastCourantNo, 0);
    addToRunTimeSelectionTable(functionObject, blastCourantNo, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::blastCourantNo::calc()
{
    if (!foundObject<volVectorField>(UName_))
    {
        return false;
    }
    const volVectorField& U = lookupObject<volVectorField>(UName_);

    tmp<volScalarField> tspeedOfSound;
    if
    (
        foundObject<fluidBlastThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        )
    )
    {
        tspeedOfSound = tmp<volScalarField>
        (
            lookupObject<fluidBlastThermo>
            (
                IOobject::groupName(physicalProperties::typeName, phaseName_)
            ).speedOfSound()
        );
    }
    else if
    (
        foundObject<fluidThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        )
    )
    {
        const fluidThermo& thermo
        (
            lookupObject<fluidThermo>
            (
                IOobject::groupName(physicalProperties::typeName, phaseName_)
            )
        );
        tspeedOfSound = sqrt(thermo.Cp()/thermo.Cv()/thermo.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Valid thermos are: "
            << mesh_.lookupClass<fluidThermo>().toc()<< endl
            << abort(FatalError);
    }

    surfaceScalarField amaxSf
    (
        fvc::interpolate(tspeedOfSound)*mesh_.magSf()
    );
    // Remove wave speed from wedge boundaries
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<wedgePolyPatch>(mesh_.boundaryMesh()[patchi]))
        {
            amaxSf.boundaryFieldRef()[patchi] = Zero;
        }
    }
    amaxSf += mag(fvc::relative(fvc::flux(U), U));

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    tmp<volScalarField> tCo
    (
        volScalarField::New
        (
            resultName_,
            mesh_,
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tCo->primitiveFieldRef() =
        0.5*mesh_.time().deltaT()
       *fvc::surfaceSum(amaxSf)()()
       /mesh_.Vsc();



    return store(resultName_, tCo);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blastCourantNo::blastCourantNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression
    (
        name,
        runTime,
        dict,
        IOobject::groupName("Ma", dict.lookupOrDefault("phaseName", word::null)),
        IOobject::groupName("U", dict.lookupOrDefault("phaseName", word::null))
    ),
    phaseName_(dict.lookupOrDefault("phaseName", word::null)),
    UName_(dict.lookupOrDefault("UName", IOobject::groupName("U", phaseName_)))
{
    if (!dict.lookupOrDefault("executeAtStart", false))
    {
        executeAtStart_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::blastCourantNo::~blastCourantNo()
{}


// ************************************************************************* //
