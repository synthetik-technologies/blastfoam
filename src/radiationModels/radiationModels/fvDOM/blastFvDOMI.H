/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

inline const Foam::radiationModels::blastRadiativeIntensityRay&
Foam::radiationModels::blastFvDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiationModels::blastFvDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiationModels::blastFvDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiationModels::blastFvDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiationModels::blastFvDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiationModels::blastFvDOM::nLambda() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::radiationModels::blastFvDOM::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::radiationModels::blastFvDOM::aLambda
(
    const label lambdaI
) const
{
    return aLambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiationModels::blastFvDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::radiationModels::blastFvDOM::qr() const
{
    return qr_;
}

inline const Foam::volScalarField& Foam::radiationModels::blastFvDOM::qin() const
{
    return qin_;
}


inline const Foam::volScalarField& Foam::radiationModels::blastFvDOM::qem() const
{
    return qem_;
}


inline const Foam::radiationModels::blastBlackBodyEmission&
Foam::radiationModels::blastFvDOM::blackBody() const
{
    return blackBody_;
}


inline Foam::scalar Foam::radiationModels::blastFvDOM::omegaMax() const
{
    return omegaMax_;
}


// ************************************************************************* //
