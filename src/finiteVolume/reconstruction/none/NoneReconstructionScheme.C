/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
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

#include "NoneReconstructionScheme.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::NoneReconstructionScheme<Type>::NoneReconstructionScheme
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    Istream& is,
    const bool overwrite
)
:
    ReconstructionScheme<Type>(phi, is, overwrite),
    tokens_(),
    is_(is.name(), tokens_)
{
    while (is.good())
    {
        tokens_.append(token(is));
    }
    is_.tokenList::operator=(tokens_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::NoneReconstructionScheme<Type>::~NoneReconstructionScheme()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::NoneReconstructionScheme<Type>::interpolateOwn() const
{
    if (!interp_.valid())
    {
        is_.rewind();
        interp_ = surfaceInterpolationScheme<Type>::New
        (
            this->phi_.mesh(),
            is_
        );
    }
    return GeometricField<Type, fvsPatchField, surfaceMesh>::New
    (
        this->ownName(),
        interp_().interpolate(this->phi_)
    );
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::NoneReconstructionScheme<Type>::interpolateNei() const
{
    if (!interp_.valid())
    {
        is_.rewind();
        interp_ = surfaceInterpolationScheme<Type>::New
        (
            this->phi_.mesh(),
            is_
        );
    }
    return GeometricField<Type, fvsPatchField, surfaceMesh>::New
    (
        this->neiName(),
        interp_().interpolate(this->phi_)
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::NoneReconstructionScheme<Type>::interpolate
(
    const surfaceScalarField& faceFlux
) const
{
    const dictionary& schemesDict =
        this->phi_.mesh().schemes().subDict("interpolationSchemes");
    return GeometricField<Type, fvsPatchField, surfaceMesh>::New
    (
        this->neiName(),
        fvc::interpolate
        (
            this->phi_,
            faceFlux,
            schemesDict.lookup("upwindInterpolate(" + this->phi_.name() +")")
        )
    );
}

// ************************************************************************* //
