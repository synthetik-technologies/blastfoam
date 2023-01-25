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
#include "hashedWordList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentumStabilisation, 0);

template<>
const char* NamedEnum<momentumStabilisation::Method, 4>::names[] =
{
    "none",
    "RhieChow",
    "Laplacian",
    "JST"
};

const NamedEnum<momentumStabilisation::Method, 4>
    momentumStabilisation::stabilisationMethods;

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentumStabilisation::momentumStabilisation
(
    const dictionary& dict
)
:
    dict_(dict.optionalSubDict("stabilisation"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentumStabilisation::setMethods(const Method defaultMethod) const
{
    setMethods(defaultMethod, dict_);
}

void Foam::momentumStabilisation::setMethods
(
    const Method defaultMethod,
    const dictionary& dict
) const
{
    const hashedWordList methods
    (
        dict.found("methods")
      ? dict.lookup<wordList>("methods")
      : wordList
        (
            1,
            dict.lookupOrDefault<word>("method", stabilisationMethods[defaultMethod])
        )
    );
    forAll(methods, i)
    {
        stabilisationMethods[methods[i]];
    }

    methods_.clear();

    // Calculate stabilisation term
    if (methods.found(stabilisationMethods[RHIE_CHOW]))
    {
        methods_.insert
        (
            RHIE_CHOW,
            dict.lookupOrDefault
            (
                word(stabilisationMethods[RHIE_CHOW]) + "ScaleFactor",
                0.1
            )
        );
    }
    if (methods.found(stabilisationMethods[LAPLACIAN]))
    {
        methods_.insert
        (
            LAPLACIAN,
            dict.lookupOrDefault
            (
                word(stabilisationMethods[LAPLACIAN]) + "ScaleFactor",
                0.01
            )
        );
    }
    if (methods.found(stabilisationMethods[JST]))
    {
        methods_.insert
        (
            JST,
            dict.lookupOrDefault
            (
                word(stabilisationMethods[JST]) + "ScaleFactor",
                0.001
            )
        );
    }
    if (methods.found(stabilisationMethods[NONE]))
    {
        methods_.insert
        (
            NONE,
            0.0
        );
    }
    Info<< "stabilisationMethods:" << incrIndent << endl;
    forAllConstIter(Map<scalar>, methods_, iter)
    {
        Method method = static_cast<Method>(iter.key());
        Info<< indent << stabilisationMethods[method] << ": " << iter() << endl;
    }
}

Foam::tmp<Foam::volVectorField> Foam::momentumStabilisation::stabilisation
(
    const volVectorField& vf,
    const volTensorField& gradVf,
    const volScalarField& gamma
) const
{
    // Lookup method
    if (!methods_.size())
    {
        setMethods(RHIE_CHOW);
    }

    tmp<volVectorField> tresult
    (
        volVectorField::New
        (
           word(type() + "Field"),
           vf.mesh(),
            dimensionedVector
            (
                "zero",
                gamma.dimensions()*gradVf.dimensions()/dimLength,
                vector::zero
            )
        )
    );

    // Calculate stabilisation term
    if (methods_.found(RHIE_CHOW))
    {
        tresult.ref() += RhieChow(methods_[RHIE_CHOW], vf, gradVf, gamma);
    }
    if (methods_.found(LAPLACIAN))
    {
        tresult.ref() += Laplacian(methods_[LAPLACIAN], vf, gamma);
    }
    if (methods_.found(JST))
    {
        tresult.ref() += JamesonSchmidtTurkel(methods_[JST], vf, gamma);
    }

    return tresult;
}


Foam::tmp<Foam::volVectorField> Foam::momentumStabilisation::stabilisation
(
    const volVectorField& vf,
    const surfaceScalarField& gamma
) const
{
    // Lookup method
    if (!methods_.size())
    {
        setMethods(JST);
    }

    tmp<volVectorField> tresult
    (
        volVectorField::New
        (
           word(type() + "Field"),
           vf.mesh(),
            dimensionedVector
            (
                "zero",
                gamma.dimensions()*vf.dimensions()/sqr(dimLength),
                vector::zero
            )
        )
    );

    // Calculate stabilisation term
    if (methods_.found(RHIE_CHOW))
    {
        FatalErrorInFunction
            << "Not implemented for velocity" << endl
            << abort(FatalError);
        // tmp<volTensorField> gradVf(fvc::grad(vf));
        // tresult.ref() += RhieChow(methods_[RHIE_CHOW], vf, gradVf(), gamma);
    }
    if (methods_.found(LAPLACIAN))
    {
        tresult.ref() += Laplacian(methods_[LAPLACIAN], vf, gamma);
    }
    if (methods_.found(JST))
    {
        tresult.ref() += JamesonSchmidtTurkel(methods_[JST], vf, gamma);
    }

    return tresult;
}


Foam::scalar Foam::momentumStabilisation::energy
(
    const volVectorField& vf,
    const volTensorField& gradVf,
    const volScalarField& gamma,
    const volTensorField& gradDD
) const
{
    NotImplemented;
    // if (!inUse())
    {
        return 0.0;
    }
/*
    tensorField smoothing(vf.mesh().nCells(), Zero);

    // Calculate stabilisation term
    if (methods_.found(RHIE_CHOW))
    {
        FatalErrorInFunction
            << "Not implemented for velocity" << endl
            << abort(FatalError);
        // smoothing += RhieChowEnergy(methods_[RHIE_CHOW], vf, gradVf, gamma);
    }
    if (methods_.found(LAPLACIAN))
    {
        smoothing += LaplacianEnergy(methods_[LAPLACIAN], vf, gamma);
    }
    if (methods_.found(JST))
    {
        smoothing += JamesonSchmidtTurkelEnergy(methods_[JST], vf, gamma);
    }

    return gSum(smoothing && gradDD.primitiveField()*vf.mesh().V());*/
}


Foam::scalar Foam::momentumStabilisation::energy
(
    const volVectorField& vf,
    const surfaceScalarField& gamma,
    const volTensorField& gradDD
) const
{
    if (!inUse())
    {
        return 0.0;
    }

    tensorField smoothing(vf.mesh().nCells(), Zero);

    // Calculate stabilisation term
    if (methods_.found(RHIE_CHOW))
    {
        // WarningInFunction
            // << "Not implemented for velocity" << endl;
        // smoothing += RhieChowEnergy(methods_[RHIE_CHOW], vf, gradVf, gamma);
    }
    if (methods_.found(LAPLACIAN))
    {
        smoothing += LaplacianEnergy(methods_[LAPLACIAN], vf, gamma);
    }
    if (methods_.found(JST))
    {
        smoothing += JamesonSchmidtTurkelEnergy(methods_[JST], vf, gamma);
    }

    return gSum(smoothing && gradDD.primitiveField()*vf.mesh().V());
}


template<>
Foam::tmp<Foam::volVectorField> Foam::momentumStabilisation::RhieChow
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

// ************************************************************************* //
