/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "sonicOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "thermodynamicConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sonicOutletVelocityFvPatchVectorField::
sonicOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    thermoBasePatchField(this->patch()),
    pRef_(Foam::constant::thermodynamic::Pstd)
{}


Foam::sonicOutletVelocityFvPatchVectorField::
sonicOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    thermoBasePatchField(p, dict),
    pRef_
    (
        dict.lookupOrDefault("pRef", Foam::constant::thermodynamic::Pstd)
    )
{}


Foam::sonicOutletVelocityFvPatchVectorField::
sonicOutletVelocityFvPatchVectorField
(
    const sonicOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    thermoBasePatchField(ptf),
    pRef_(ptf.pRef_)
{}


Foam::sonicOutletVelocityFvPatchVectorField::
sonicOutletVelocityFvPatchVectorField
(
    const sonicOutletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    thermoBasePatchField(ptf),
    pRef_(ptf.pRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sonicOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField gamma(this->gamma());
    const scalarField& p =
        this->thermo().p().boundaryField()[this->patch().index()];
//     const scalar T =
//         gSum
//         (
//             this->thermo().T().boundaryField()[this->patch().index()]
//            *magSf
//         )/A;
//     const scalar R =
//         gSum
//         (
//             Foam::constant::thermodynamic::RR
//            /this->thermo().W(this->patch().index())
//            *magSf
//         )/A;

    scalarField Ma(this->size(), 0.0);
    forAll(*this, i)
    {
        const scalar g = gamma[i];
        scalar pExit = p[i]*pow(2.0/(g + 1.0), g/(g - 1.0));
//         scalar Texit = T;
//         scalar G = 0.0;
        if (pExit >= pRef_)
        {
            Ma = 1.0;
//             Texit = T*2.0/(g + 1.0);
//             G =
//                 p
//                *sqrt
//                 (
//                     pow(2.0/(g + 1.0), (g + 1.0)/(g - 1.0))
//                    *(g/R/T)
//                 );
        }
        else if (p[i] > pRef_)
        {
            pExit = pRef_;
            Ma[i] =
                sqrt
                (
                    max
                    (
                        (pow(p[i]/pRef_, (g - 1.0)/g) - 1.0)
                       *2.0/(g - 1.0),
                        0.0
                    )
                );
//             Texit = T/(1.0 + sqr(Ma)*(gamma - 1.0)/2.0);
//             G =
//                 p
//                *sqrt(gamma/R/T)
//                *Ma
//                /pow
//                 (
//                     1.0 + (gamma - 1.0)/2.0*sqr(Ma),
//                     (gamma + 1.0)/(2.0*(gamma - 1.0))
//                 );
        }
    }
    const scalarField& magSf = this->patch().magSf();
    const scalar A = gSum(magSf);
    tmp<scalarField> tc(this->speedOfSound());
    Info<< this->patch().name() << ":" << nl
        << "    Mach number: " << gSum(magSf*Ma)/A << nl
        << "    Speed of sound: " << gSum(tc()*magSf)/A << nl
        << "    Exit velocity: " << gSum(Ma*tc()*magSf)/A << nl << endl;
    vectorField normal(this->patch().nf());
    vectorField UI(this->patchInternalField());
    scalarField Un(UI & normal);
    vectorField Ut(UI - Un*normal);
    vectorField::operator=(Ut + Ma*tc*normal);
    fixedValueFvPatchVectorField::updateCoeffs();
}



void Foam::sonicOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    thermoBasePatchField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        sonicOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
