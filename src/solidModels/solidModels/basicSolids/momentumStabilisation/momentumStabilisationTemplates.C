/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "momentumStabilisation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GammaField>
Foam::tmp<Foam::volVectorField> Foam::momentumStabilisation::RhieChow
(
    const scalar scale,
    const volVectorField& vf,
    const volTensorField& gradVf,
    const GammaField& gamma
) const
{
    return
        scale
       *(
            fvc::laplacian
            (
                gamma,
                vf,
                "laplacian(D" + vf.name() +"," + vf.name() + ")"
            )
          - fvc::div(gamma*gradVf)
        );
}


template<class GammaField>
Foam::tmp<Foam::volVectorField> Foam::momentumStabilisation::JamesonSchmidtTurkel
(
    const scalar scale,
    const volVectorField& vf,
    const GammaField& gamma
) const
{
    const word scheme("laplacian(D" + vf.name() +"," + vf.name() + ")");
    return
      - scale
       *fvc::laplacian
        (
            vf.mesh().magSf(),
            fvc::laplacian(gamma, vf, scheme),
            scheme
        );
}


template<class GammaField>
Foam::tmp<Foam::volVectorField> Foam::momentumStabilisation::Laplacian
(
    const scalar scale,
    const volVectorField& vf,
    const GammaField& gamma
) const
{
    return scale*fvc::laplacian(gamma, vf);
}


template<class GammaField>
Foam::tmp<Foam::tensorField>
Foam::momentumStabilisation::RhieChowEnergy
(
    const scalar scale,
    const volVectorField& vf,
    const volTensorField& gradVf,
    const GammaField& gamma
) const
{
    return
        scale
       *(
            fvc::reconstruct
            (
                gamma*0.5
               *(
                   fvc::snGrad(vf) + fvc::snGrad(vf.oldTime())
                )*vf.mesh().magSf()
            )().primitiveField()
          - fvc::div(gamma*gradVf)().primitiveField()
        );
}


template<class GammaField>
Foam::tmp<Foam::tensorField>
Foam::momentumStabilisation::JamesonSchmidtTurkelEnergy
(
    const scalar scale,
    const volVectorField& vf,
    const GammaField& gamma
) const
{
    const word scheme("laplacian(D" + vf.name() +"," + vf.name() + ")");
    tmp<volVectorField> lapGammaVf(fvc::laplacian(gamma, vf, scheme));
    tmp<volVectorField> lapGammaVfOld
    (
        fvc::laplacian(gamma, vf.oldTime(), scheme)
    );
    return
      - scale
       *fvc::reconstruct
        (
            0.5*vf.mesh().magSf()
           *(
                fvc::snGrad(lapGammaVf)
              + fvc::snGrad(lapGammaVfOld)
            )
        )().primitiveField();
}


template<class GammaField>
Foam::tmp<Foam::tensorField>
Foam::momentumStabilisation::LaplacianEnergy
(
    const scalar scale,
    const volVectorField& vf,
    const GammaField& gamma
) const
{
    return
        scale
       *(
            fvc::reconstruct
            (
                gamma*0.5
               *(
                   fvc::snGrad(vf) + fvc::snGrad(vf.oldTime())
                )*vf.mesh().magSf()
            )().primitiveField()
        );
}

// ************************************************************************* //
