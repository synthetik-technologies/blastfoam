/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
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

Class
    Foam::compressible::correction

Description
    Implementation of the Sarkar/Zemen compressibility correction

    \verbatim
        Sarkar, S., Erlebacher, G., Hussaini, M. Y., and Kreiss, H. O. (2001).
        The Analysis and Modeling of Dilatational Terms in Compressible Turbulence.
        NASA CR 181959
    \endverbatim

    \verbatim
         Zeman, O. (1990).
        Dilatational Dissipation: The Concept and Application in Modeling
        Compressible Mixing Layers.
        Physics of Fluids A, Vol. 2, No. 2,
    \endverbatim

SourceFiles
    compressibilityCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef compressibilityCorrection_H
#define compressibilityCorrection_H

#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
                          Class compressible::correction Declaration
\*---------------------------------------------------------------------------*/

class correction
{
protected:

    // Protected data

        // Model coefficients

            dimensionedScalar XiStar_;

            dimensionedScalar Mat0_;

            dimensionedScalar rProdLim_;


public:

    // Constructors

        //- Construct from components
        correction(dictionary& dict);

        //- Disallow default bitwise copy construction
        correction(const correction&) = delete;


    //- Destructor
    virtual ~correction()
    {}

    // Member Functions

        //- Read RASProperties dictionary
        bool read(const dictionary&);

        //- Return tubulent Mach number
        static tmp<volScalarField::Internal> Mt
        (
            const volScalarField& k
        );

        //- Return tubulent Mach number squared
        static tmp<volScalarField::Internal> MtSqr
        (
            const volScalarField& k
        );

        //- Return tubulent Mach number squared divided by k
        //  for use with implicit sources
        static tmp<volScalarField::Internal> MtSqrByk
        (
            const volScalarField& k
        );

        //- Return compressible correction factor
        tmp<volScalarField::Internal> Fcorr
        (
            const volScalarField& k
        ) const;

        //- Return compressible beta star
        void correct
        (
            const volScalarField& k,
            const dimensionedScalar& betaStar0,
            tmp<volScalarField::Internal>& beta,
            tmp<volScalarField::Internal>& betaStar
        ) const;

        //- Return compressible beta star
        void correct
        (
            const volScalarField& k,
            const dimensionedScalar& beta0,
            const dimensionedScalar& betaStar0,
            tmp<volScalarField::Internal>& beta,
            tmp<volScalarField::Internal>& betaStar
        ) const;

        //- Return limited production source
        void limitG
        (
            volScalarField::Internal& G,
            const volScalarField::Internal& k,
            const volScalarField::Internal& omega,
            const volScalarField::Internal& betaStar
        ) const;

        //- Limit epsilon based on length scale
        static void limitEpsilon_L
        (
            const volScalarField& k,
            const volScalarField& y,
            volScalarField& epsilon
        );

        //- Limit omega based on length scale
        static void limitOmega_L
        (
            const volScalarField& k,
            const volScalarField& y,
            volScalarField& omega
        );
};


#define CompressibleTypeName(Name)                                             \
    TypeName("compressible::" Name);                                           \
    static word baseName() {return Name;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
