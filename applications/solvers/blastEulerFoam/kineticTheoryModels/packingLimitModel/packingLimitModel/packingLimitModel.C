/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "packingLimitModel.H"
#include "SortableList.H"
#include "constantDiameterModel.H"
#include "zeroGradientFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(packingLimitModel, 0);

    defineRunTimeSelectionTable(packingLimitModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModel::packingLimitModel
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    dict_(dict),
    kt_(kt),
    mesh_(kt.fluid().mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModel::~packingLimitModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::packingLimitModel::alphaMax() const
{
    const UPtrList<const phaseModel>& phases(kt_.packingPhases());
    tmp<volScalarField> tmpAlphaMax
    (
        new volScalarField
        (
            IOobject
            (
                "alphaMax",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "alphaMax",
                dimless,
                phases[0].alphaMax()
            ),
            wordList
            (
                mesh_.boundaryMesh().size(),
                zeroGradientFvPatchScalarField::typeName
            )
        )
    );
    volScalarField& alphaMaxField(tmpAlphaMax.ref());

    if (phases.size() == 1)
    {
        return tmpAlphaMax;
    }


    bool constantDiameters = true;
    forAll(phases, phasei)
    {
        if (!isA<diameterModels::constantDiameter>(phases[phasei].dModel()))
        {
            constantDiameters = false;
        }
    }

    // Only sort diameters in one cell to save time
    if (constantDiameters)
    {
        // Sort diameters from largest to smallest
        scalarList ds(phases.size());
        forAll(phases, phasei)
        {
            ds[phasei] = phases[phasei].d()()[0];
        }

        alphaMaxField.primitiveFieldRef() = alphaMax(0, ds);
    }
    // Sort particle diameters for every cell
    else
    {
        forAll(alphaMaxField, celli)
        {
            // Sort diameters from largest to smallest
            scalarList ds(phases.size());

            forAll(phases, phasei)
            {
                ds[phasei] = phases[phasei].celld(celli);
            }

            alphaMaxField[celli] = alphaMax(celli, ds);
        }
    }
    alphaMaxField.correctBoundaryConditions();

    return tmpAlphaMax;
}
// ************************************************************************* //
