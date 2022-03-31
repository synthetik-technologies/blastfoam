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
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentumStabilisation, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentumStabilisation::momentumStabilisation
(
    const dictionary& dict
)
:
    dict_(dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::momentumStabilisation::stabilisation
(
    const volVectorField& vf,
    const volTensorField& gradVf,
    const volScalarField& gamma
) const
{
    tmp<volVectorField> tresult
    (
        new volVectorField
        (
            IOobject
            (
                word(type() + "Field"),
                vf.mesh().time().timeName(),
                vf.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            dimensionedVector
            (
                "zero",
                gamma.dimensions()*gradVf.dimensions()/dimLength,
                vector::zero
            )
        )
    );
    volVectorField& result = tresult.ref();

    // Lookup method
    const word method = word(dict_.lookup("type"));

    // Calculate stabilisation term
    if (method == "RhieChow")
    {
        const scalar scaleFactor = readScalar(dict_.lookup("scaleFactor"));

        result = scaleFactor
       *(
           fvc::laplacian(gamma, vf, "laplacian(DD,D)") - fvc::div(gamma*gradVf)
        );
    }
    else if (method == "JamesonSchmidtTurkel")
    {
        const scalar scaleFactor = readScalar(dict_.lookup("scaleFactor"));

        result = -scaleFactor*fvc::laplacian
        (
            vf.mesh().magSf(),
            //1.0/(vf.mesh().deltaCoeffs()*vf.mesh().deltaCoeffs()),
            fvc::laplacian(gamma, vf, "JSTinner"),
            "JSTouter"
        );
    }
    else if (method == "Laplacian")
    {
        const scalar scaleFactor = readScalar(dict_.lookup("scaleFactor"));

        result = scaleFactor*fvc::laplacian(gamma, vf);
    }
    else if (method != "none")
    {
        FatalErrorIn(type() + "::stabilisation() const")
            << "Unknown method = " << method << nl
            << "Methods are: none, RhieChow, JamesonSchmidtTurkel and Laplacian"
            <<  abort(FatalError);
    }

    return tresult;
}



// ************************************************************************* //
