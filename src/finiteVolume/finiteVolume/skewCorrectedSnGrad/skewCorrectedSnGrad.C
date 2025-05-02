/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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

#include "skewCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::skewCorrectedSnGrad<Type>::fullGradCorrection
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            "snGradCorr("+vf.name()+')',
            mesh,
            vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf.ref();
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bssf =
        ssf.boundaryFieldRef();
    ssf = dimensioned<Type>(ssf.dimensions(), Zero);


    typedef typename
        outerProduct<vector, typename pTraits<Type>::cmptType>::type
        CmptGradType;

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const vectorField& Sf = mesh.Sf().primitiveField();
    const scalarField& magSf = mesh.magSf().primitiveField();

    vectorField nf(Sf/magSf);

    const vectorField& Cf = mesh.Cf().primitiveField();
    const vectorField& C = mesh.C().primitiveField();

    const scalarField& deltaCoeffs =
        mesh.deltaCoeffs().primitiveField();

    surfaceVectorField kP
    (
        IOobject
        (
            "kP",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimLength, Zero)
    );
    vectorField& kPI = kP.primitiveFieldRef();
    surfaceVectorField::Boundary& bkP = kP.boundaryFieldRef();

    surfaceVectorField kN
    (
        IOobject
        (
            "kN",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimLength, Zero)
    );
    vectorField& kNI = kN.primitiveFieldRef();
    surfaceVectorField::Boundary& bkN = kN.boundaryFieldRef();

    kPI = Cf - vectorField(C, owner);
    kPI -= Sf*(Sf & kPI)/sqr(magSf);

    kNI = Cf - vectorField(C, neighbour);
    kNI -= Sf*(Sf & kNI)/sqr(magSf);

    forAll(bkP, patchI)
    {
        if (bkP[patchI].coupled())
        {
            const fvPatch& patch = mesh.boundary()[patchI];

            bkP[patchI] =
                patch.Cf() - patch.Cn()
              - patch.Sf()*(patch.Sf() & bkP[patchI])/sqr(patch.magSf());

            bkN[patchI] =
                patch.Cf() - (patch.Cn() + patch.delta())
              - patch.Sf()*(patch.Sf() & bkN[patchI])/sqr(patch.magSf());
        }
    }

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        GeometricField<CmptGradType, fvPatchField, volMesh> cmptGrad
        (
            gradScheme<typename pTraits<Type>::cmptType>::New
            (
                mesh,
                mesh.schemes().grad("grad(" + vf.name() + ')')
            )()
           .grad(vf.component(cmpt))
        );

        const Field<CmptGradType>& cmptGradI = cmptGrad.primitiveField();

        // Skewness and non-rothogonal correction
        {
            ssf.primitiveFieldRef().replace
            (
                cmpt,
                (
                    (kNI & Field<CmptGradType>(cmptGradI, neighbour))
                  - (kPI & Field<CmptGradType>(cmptGradI, owner))
                )
               *deltaCoeffs
            );
        }

        forAll(bssf, patchI)
        {
            if (bssf[patchI].coupled())
            {
                bssf[patchI].replace
                (
                    cmpt,
                    (
                        (
                            bkN[patchI]
                          & cmptGrad.boundaryField()[patchI].patchNeighbourField()
                        )
                      - (
                            bkP[patchI]
                          & cmptGrad.boundaryField()[patchI].patchInternalField()
                        )
                    )
                   *mesh.deltaCoeffs().boundaryField()[patchI]
                );
            }
        }
    }

    // // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    // tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf =
    //     linear<typename outerProduct<vector, Type>::type>(mesh).dotInterpolate
    //     (
    //         mesh.nonOrthCorrectionVectors(),
    //         gradScheme<Type>::New
    //         (
    //             mesh,
    //             mesh.gradScheme("grad(" + vf.name() + ')')
    //         )().grad(vf, "grad(" + vf.name() + ')')
    //     );
    //
    // tssf.ref().rename("snGradCorr(" + vf.name() + ')');

    return tssf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::skewCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            "snGradCorr(" + vf.name() + ')',
            mesh,
            vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf.ref();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        ssf.replace
        (
            cmpt,
            skewCorrectedSnGrad<typename pTraits<Type>::cmptType>
            (
                mesh
            ).fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}



// ************************************************************************* //
