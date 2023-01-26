/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "dynamicMultiAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

dynamicMultiAlphaContactAngleFvPatchScalarField::interfaceThetaProps::interfaceThetaProps
(
    Istream& is
)
:
    theta0_(readScalar(is)),
    uTheta_(readScalar(is)),
    thetaA_(readScalar(is)),
    thetaR_(readScalar(is))
{}


Istream& operator>>
(
    Istream& is,
    dynamicMultiAlphaContactAngleFvPatchScalarField::interfaceThetaProps& tp
)
{
    is >> tp.theta0_ >> tp.uTheta_ >> tp.thetaA_ >> tp.thetaR_;
    return is;
}


Ostream& operator<<
(
    Ostream& os,
    const dynamicMultiAlphaContactAngleFvPatchScalarField::interfaceThetaProps& tp
)
{
    os  << tp.theta0_ << token::SPACE
        << tp.uTheta_ << token::SPACE
        << tp.thetaA_ << token::SPACE
        << tp.thetaR_;

    return os;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicMultiAlphaContactAngleFvPatchScalarField::
dynamicMultiAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF)
{}


dynamicMultiAlphaContactAngleFvPatchScalarField::
dynamicMultiAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    thetaProps_(dict.lookup("thetaProperties"))
{
    evaluate();
}


dynamicMultiAlphaContactAngleFvPatchScalarField::
dynamicMultiAlphaContactAngleFvPatchScalarField
(
    const dynamicMultiAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    thetaProps_(gcpsf.thetaProps_)
{}


dynamicMultiAlphaContactAngleFvPatchScalarField::
dynamicMultiAlphaContactAngleFvPatchScalarField
(
    const dynamicMultiAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, iF),
    thetaProps_(gcpsf.thetaProps_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynamicMultiAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleFvPatchScalarField::write(os);
    if (thetaProps_.size())
    {
        writeEntry(os, "thetaProperties", thetaProps_);
    }
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    dynamicMultiAlphaContactAngleFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
