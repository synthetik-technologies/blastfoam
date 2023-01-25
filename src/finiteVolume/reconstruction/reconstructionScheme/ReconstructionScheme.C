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

#include "ReconstructionScheme.H"
#include "StandardReconstructionScheme.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::ownName() const
{
    return this->ownName(phi_.name());
}


template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::neiName() const
{
    return this->neiName(phi_.name());
}


template<class Type>
void
Foam::ReconstructionScheme<Type>::interpolateOwnNei
(
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tphiOwn,
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tphiNei
) const
{
    tphiOwn.clear();
    if
    (
        !overwrite_
     && phi_.mesh().template foundObject
        <
            GeometricField<Type, fvsPatchField, surfaceMesh>
        >(ownName())
    )
    {
        DebugInfo << "Reading " << ownName() << " from cache" << endl;
        tphiOwn = tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        (
            phi_.mesh().template lookupObject
            <
                GeometricField<Type, fvsPatchField, surfaceMesh>
            >(ownName())
        );
    }
    else
    {
        DebugInfo << "Recomputing " << ownName() << endl;
        tphiOwn = interpolateOwn();
    }

    tphiNei.clear();
    if
    (
        !overwrite_
     && phi_.mesh().template foundObject
        <
            GeometricField<Type, fvsPatchField, surfaceMesh>
        >(neiName())
    )
    {
        DebugInfo << "Reading " << neiName() << " from cache" << endl;
        tphiNei = tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        (
            phi_.mesh().template lookupObject
            <
                GeometricField<Type, fvsPatchField, surfaceMesh>
            >(neiName())
        );
    }
    else
    {
        DebugInfo << "Recomputing " << neiName() << endl;
        tphiNei = interpolateNei();
    }
}


template<class Type>
Foam::autoPtr<Foam::ReconstructionScheme<Type>>
Foam::ReconstructionScheme<Type>::New
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const word& fieldName,
    const bool overwrite
)
{
    return New(phi, fieldName, word::null, overwrite);
}


template<class Type>
Foam::autoPtr<Foam::ReconstructionScheme<Type>>
Foam::ReconstructionScheme<Type>::New
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const word& fieldName,
    const word& phaseName,
    const bool overwrite
)
{
    const word schemeKey(scheme(fieldName, phaseName, phi.mesh(), debug, overwrite));
    Istream& is(phi.mesh().interpolationScheme(schemeKey));
    word order(is);
    word scheme(order);
    typedef surfaceInterpolationScheme<Type> sISType;
    typename sISType::MeshFluxConstructorTable::iterator iter =
        sISType::MeshFluxConstructorTablePtr_->find(scheme);

    // Standard OpenFOAM interpolation
    if (iter != sISType::MeshFluxConstructorTablePtr_->end())
    {
        return autoPtr<ReconstructionScheme<Type>>
        (
            new StandardReconstructionScheme<Type>
            (
                phi,
                IStringStream(schemeKey)(),
                overwrite
            )
        );
    }

    if (is.good() && scheme != "none" && scheme != "upwindMUSCL" && scheme != "THINC")
    {
        token t(is);
        if (t.isWord())
        {
            scheme = scheme + '<' + t.wordToken() + '>';
        }
        else
        {
            is.putBack(t);
        }
    }

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(scheme);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Unknown discretisation scheme "
            << scheme << " for " << fieldName << nl << nl
            << "Valid MUSCL schemes are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl << nl
            << "Valid OpenFOAM schemes are:" << nl
            << sISType::MeshFluxConstructorTablePtr_->sortedToc()
            << endl
            << abort(FatalIOError);
    }

    return cstrIter()(phi, is, overwrite);
}\

// ************************************************************************* //
