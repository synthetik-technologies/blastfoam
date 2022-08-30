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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "analyticalPlateHoleTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mechanicalModel.H"
#include "volFields.H"
#include "fvc.H"
#include "fixedValueFvPatchFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

symmTensor analyticalPlateHoleTractionFvPatchVectorField::plateHoleSolution
(
    const vector& C
)
{
    tensor sigma = tensor::zero;

    // Calculate radial coordinate
    scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));

    // Calculate circumferential coordinate
    scalar theta = Foam::atan2(C.y(), C.x());

    coordinateSystem cs("polarCS", C, vector(0, 0, 1), C/mag(C));

    sigma.xx() =
        T_*(1 - sqr(holeR_)/sqr(r))/2
      + T_
       *(1 + 3*pow(holeR_,4)/pow(r,4) - 4*sqr(holeR_)/sqr(r))*::cos(2*theta)/2;

    sigma.xy() =
      - T_
       *(1 - 3*pow(holeR_,4)/pow(r,4) + 2*sqr(holeR_)/sqr(r))*::sin(2*theta)/2;

    sigma.yx() = sigma.xy();

    sigma.yy() =
        T_*(1 + sqr(holeR_)/sqr(r))/2
      - T_*(1 + 3*pow(holeR_,4)/pow(r,4))*::cos(2*theta)/2;


    // Transformation to Cartesian coordinate system
    sigma = ((cs.R().R() & sigma) & cs.R().R().T());

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();

    return S;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    T_(0.0),
    holeR_(0.0)
{}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    T_(readScalar(dict.lookup("farFieldTractionX"))),
    holeR_(readScalar(dict.lookup("holeRadius")))
{}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(stpvf, p, iF, mapper),
    T_(stpvf.T_),
    holeR_(stpvf.holeR_)
{}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(stpvf, iF),
    T_(stpvf.T_),
    holeR_(stpvf.holeR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void analyticalPlateHoleTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void analyticalPlateHoleTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void analyticalPlateHoleTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Patch unit normals
    vectorField n(patch().nf());

    // Patch face centres
    const vectorField& Cf = patch().Cf();

    // Set the patch traction

    vectorField& trac = traction();

    forAll(traction(), faceI)
    {
        vector curC(Cf[faceI].x(), Cf[faceI].y(), 0);
        vector curN = n[faceI];

        if (patch().name() == "hole")
        {
            curC /= mag(curC);
            curC *= holeR_;

            curN = -curC/mag(curC);
        }

        trac[faceI] = (n[faceI] & plateHoleSolution(curC));
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}


// Write
void analyticalPlateHoleTractionFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("farFieldTractionX")
        << T_ << token::END_STATEMENT << nl;

    os.writeKeyword("holeRadius")
        << holeR_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    analyticalPlateHoleTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
