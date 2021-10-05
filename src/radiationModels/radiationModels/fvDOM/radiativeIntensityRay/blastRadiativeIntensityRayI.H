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

inline const Foam::volScalarField&
Foam::radiationModels::blastRadiativeIntensityRay::I() const
{
    return I_;
}


inline const Foam::volScalarField&
Foam::radiationModels::blastRadiativeIntensityRay::qr() const
{
    return qr_;
}


inline Foam::volScalarField&
Foam::radiationModels::blastRadiativeIntensityRay::qr()
{
    return qr_;
}

inline const Foam::volScalarField&
Foam::radiationModels::blastRadiativeIntensityRay::qin() const
{
    return qin_;
}


inline Foam::volScalarField&
Foam::radiationModels::blastRadiativeIntensityRay::qin()
{
    return qin_;
}


inline const Foam::volScalarField& Foam::radiationModels::
blastRadiativeIntensityRay::qem() const
{
    return qem_;
}


inline Foam::volScalarField&
Foam::radiationModels::blastRadiativeIntensityRay::qem()
{
    return qem_;
}


inline const Foam::vector&
Foam::radiationModels::blastRadiativeIntensityRay::d() const
{
    return d_;
}


inline const Foam::vector&
Foam::radiationModels::blastRadiativeIntensityRay::dAve() const
{
    return dAve_;
}


inline Foam::scalar
Foam::radiationModels::blastRadiativeIntensityRay::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiationModels::blastRadiativeIntensityRay::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiationModels::blastRadiativeIntensityRay::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiationModels::blastRadiativeIntensityRay::omega() const
{
    return omega_;
}


inline const Foam::volScalarField&
Foam::radiationModels::blastRadiativeIntensityRay::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
