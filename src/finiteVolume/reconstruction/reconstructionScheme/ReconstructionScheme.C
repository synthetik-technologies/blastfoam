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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::scheme(const word& name)
{
    return "reconstruct(" + name + ")";
}


template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::scheme
(
    const word& name,
    const fvMesh& mesh,
    const bool fail
)
{
    return scheme(name, word::null, mesh, fail);
}


template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::scheme
(
    const word& baseName,
    const word& phaseName,
    const fvMesh& mesh,
    const bool fail
)
{
    const word name(IOobject::groupName(name, phaseName));
    word baseScheme(scheme(baseName));
    word nameScheme(scheme(name));

    if (mesh.schemesDict().subDict("interpolationSchemes").found(nameScheme))
    {
        return nameScheme;
    }
    else if (mesh.schemesDict().subDict("interpolationSchemes").found(baseScheme))
    {
        return baseScheme;
    }
    else if (fail)
    {
        FatalErrorInFunction
            << "Riemann fluxes are used, but no limiter is " << nl
            << "specified for " << name << "." << nl
            << "Please specify " << string(nameScheme)
            << " or " << string(baseScheme) << endl
            << "This may result in unstable solutions." << endl;
    }
    else
    {
        WarningInFunction
            << "Riemann fluxes are used, but no limiter is " << nl
            << "specified for " << name << "." << nl
            << "This may result in unstable solutions." << nl
            << "Please specify " << string(nameScheme)
            << " or " << string(baseScheme) << endl;

    }
    return nameScheme;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::ownName() const
{
    return ownName(phi_.name());
}


template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::ownName(const word& name) const
{
    return
        IOobject::groupName
        (
            IOobject::member(name) + "Own",
            IOobject::group(name)
        );
}


template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::neiName() const
{
    return neiName(phi_.name());
}

template<class Type>
Foam::word Foam::ReconstructionScheme<Type>::neiName(const word& name) const
{
    return
        IOobject::groupName
        (
            IOobject::member(name) + "Nei",
            IOobject::group(name)
        );
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
        DebugInfo << "Recomputing " << neiName() <<" "<<overwrite_<< endl;
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
    word name
    (
        "reconstruct("
      + fieldName
      + ")"
    );
    word namePhase
    (
        "reconstruct("
      + IOobject::groupName(fieldName, phaseName)
      + ")"
    );

    if
    (
        phi.mesh().schemesDict().subDict
        (
            "interpolationSchemes"
        ).found(namePhase)
    )
    {
        name = namePhase;
    }
    else if
    (
        !phi.mesh().schemesDict().subDict
        (
            "interpolationSchemes"
        ).found(name)
    )
    {
        // Default to lookup of density scheme
        if
        (
            !phi.mesh().schemesDict().subDict
            (
                "interpolationSchemes"
            ).found("reconstruct(rho)")
        )
        {
            WarningInFunction
                << "Riemann fluxes are used, but no limiter is " << nl
                << "specified for " << name << "." << nl
                << "This may result in unstable solutions." << endl;
        }
        name = "reconstruct(rho)";
    }

    Istream& is(phi.mesh().interpolationScheme(name));
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
                IStringStream(name)(),
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
            << scheme << " for " << name << nl << nl
            << "Valid MUSCL schemes are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc() << nl << nl
            << "Valid OpenFOAM schemes are:" << endl
            << sISType::MeshFluxConstructorTablePtr_->sortedToc()
            << abort(FatalIOError);
    }

    return cstrIter()(phi, is, overwrite);
}


// template<class Type>
// Foam::autoPtr<Foam::ReconstructionScheme<Type>>
// Foam::ReconstructionScheme<Type>::New
// (
//     const GeometricField<Type, fvPatchField, volMesh>& phi,
//     const word& fieldName,
//     const word& phaseName
// )
// {
//     word name
//     (
//         "reconstruct("
//       + fieldName
//       + ")"
//     );
//     word namePhase
//     (
//         "reconstruct("
//       + IOobject::groupName(fieldName, phaseName)
//       + ")"
//     );
//
//     if
//     (
//         phi.mesh().schemesDict().subDict
//         (
//             "interpolationSchemes"
//         ).found(namePhase)
//     )
//     {
//         name = namePhase;
//     }
//     else if
//     (
//         !phi.mesh().schemesDict().subDict
//         (
//             "interpolationSchemes"
//         ).found(name)
//     )
//     {
//         //- Default to lookup of density scheme
//         if
//         (
//             debug
//          || !phi.mesh().schemesDict().subDict
//             (
//                 "interpolationSchemes"
//             ).found("reconstruct(rho)")
//         )
//         {
//             WarningInFunction
//                 << "Riemann fluxes are used, but no limiter is " << nl
//                 << "specified for " << name << "." << nl
//                 << "This may result in unstable solutions." << endl;
//         }
//         name = "reconstruct(rho)";
//     }
//
//     Istream& is(phi.mesh().interpolationScheme(name));
//     word order(is);
//
//     if (debug)
//     {
//         Info<< "selecting " << order << " interpolation scheme "
//             << "for " << phi.name() << endl;
//     }
//
//
//     // No upwinding scheme
//     if (order == "none")
//     {
//         return autoPtr<ReconstructionScheme<Type>>
//         (
//             new noneMUSCLReconstructionScheme<Type>(phi, is)
//         );
//     }
//
//     // Upwind scheme
//     if (order == "upwindMUSCL")
//     {
//         return autoPtr<ReconstructionScheme<Type>>
//         (
//             new upwindMUSCLReconstructionScheme<Type>(phi, is)
//         );
//     }
//
//     // Linear MUSCL
//     if (order == "linearMUSCL")
//     {
//         word limiterName(is);
//         typename linearMeshConstructorTable::iterator cstrIter =
//             linearMeshConstructorTablePtr_->find(limiterName);
//
//         if (cstrIter == linearMeshConstructorTablePtr_->end())
//         {
//             FatalErrorInFunction
//                 << "Unknown linear MUSCL limiter type "
//                 << limiterName << " for " << phi.name() << endl << endl
//                 << "Valid linear MUSCL limiters types are : " << endl
//                 << linearMeshConstructorTablePtr_->sortedToc()
//                 << exit(FatalError);
//         }
//
//         return cstrIter()(phi, is);
//     }
//
//     // Quadratic MUSCL
//     if (order == "quadraticMUSCL")
//     {
//         word limiterName(is);
//         typename quadraticMeshConstructorTable::iterator cstrIter =
//             quadraticMeshConstructorTablePtr_->find(limiterName);
//
//         if (cstrIter == quadraticMeshConstructorTablePtr_->end())
//         {
//             FatalErrorInFunction
//                 << "Unknown quadratic MUSCL limiter type "
//                 << limiterName << " for " << phi.name() << endl << endl
//                 << "Valid quadratic MUSCL limiters are : " << endl
//                 << quadraticMeshConstructorTablePtr_->sortedToc()
//                 << exit(FatalError);
//         }
//
//         return cstrIter()(phi, is);
//     }
//
//     if (debug)
//     {
//         Info<< "No MUSCL scheme named " << name << ". Using standard " << nl
//             << "interpolation schemes instead." << endl;
//     }
//
//     // Standard OpenFOAM interpolation
//     IStringStream nameIS(name);
//     return autoPtr<ReconstructionScheme<Type>>
//     (
//         new standardMUSCLReconstructionScheme<Type>(phi, nameIS)
//     );
// }

// ************************************************************************* //
