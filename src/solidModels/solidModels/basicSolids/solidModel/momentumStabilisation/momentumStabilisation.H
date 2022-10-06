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

Class
    momentumStabilisation

Description
    Run-time selectable stabilisation (diffusion) term for the solid mechanics
    momentum equation.

Author
    Philip Cardiff, UCD.

SourceFiles
    momentumStabilisation.C

\*---------------------------------------------------------------------------*/

#ifndef momentumStabilisation_H
#define momentumStabilisation_H

#include "dictionary.H"
#include "volFields.H"
#include "fvc.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class momentumStabilisation Declaration
\*---------------------------------------------------------------------------*/

class momentumStabilisation
{
public:
    enum Method
    {
        NONE,
        RHIE_CHOW,
        LAPLACIAN,
        JST,
        INVALID
    };
    static const NamedEnum<Method, 4> stabilisationMethods;

private:
    // Private data

        //- Copy of settings dict
        const dictionary& dict_;

        mutable Map<scalar> methods_;


    // Private Member Functions

        template<class GammaField>
        tmp<volVectorField> RhieChow
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
        template<>
        tmp<volVectorField> RhieChow
        (
            const scalar scale,
            const volVectorField& vf,
            const volTensorField& gradVf,
            const surfaceScalarField& gamma
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
                  - fvc::div
                    (
                        gamma
                       *(
                            fvc::interpolate(gradVf) & vf.mesh().Sf()
                        )
                    )
                );
        }

        template<class GammaField>
        tmp<volVectorField> JamesonSchmidtTurkel
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
                    //1.0/(vf.mesh().deltaCoeffs()*vf.mesh().deltaCoeffs()),
                    fvc::laplacian(gamma, vf, scheme),
                    scheme
            );
        }
        template<class GammaField>
        tmp<volVectorField> Laplacian
        (
            const scalar scale,
            const volVectorField& vf,
            const GammaField& gamma
        ) const
        {
            return scale*fvc::laplacian(gamma, vf);
        }

        //- Disallow default bitwise copy construct
        momentumStabilisation(const momentumStabilisation&);

        //- Disallow default bitwise assignment
        void operator=(const momentumStabilisation&);

public:

    //- Runtime type information
    TypeName("momentumStabilisation");


    // Constructors

        //- Construct from dictionary
        momentumStabilisation
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~momentumStabilisation()
        {}


    // Member Functions

        //- Set the stabilisation methods
        void setMethods(const Method defaultMethod) const;
        void setMethods(const Method defaultMethod, const dictionary& dict) const;

        // Evaluation

            //- Stabilisation term
            tmp<volVectorField> stabilisation
            (
                const volVectorField& vf,
                const volTensorField& gradVf,
                const volScalarField& gamma
            ) const;

            //- Stabilisation term
            tmp<volVectorField> stabilisation
            (
                const volVectorField& vf,
                const volTensorField& gradVf,
                const surfaceScalarField& gamma
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
