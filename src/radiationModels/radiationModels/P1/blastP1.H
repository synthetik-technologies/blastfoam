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

Class
    Foam::radiationModels::blastP1

Description
    Works well for combustion applications where optical thickness, tau is
    large, i.e. tau = a*L > 3 (L = distance between objects)

    Assumes
     - all surfaces are diffuse
     - tends to over predict radiative fluxes from sources/sinks
       *** SOURCES NOT CURRENTLY INCLUDED ***

SourceFiles
    blastP1.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_blastP1_H
#define radiation_blastP1_H

#include "blastRadiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                           Class blastP1 Declaration
\*---------------------------------------------------------------------------*/

class blastP1
:
    public blastRadiationModel
{
    // Private Data

        //- Incident radiation / [W/m^2]
        volScalarField G_;

        //- Total radiative heat flux [W/m^2]
        volScalarField qr_;

        //- Absorption coefficient
        volScalarField a_;

        //- Emission coefficient
        volScalarField e_;

        //- Emission contribution
        volScalarField E_;


public:

    //- Runtime type information
    TypeName("P1");


    // Constructors

        //- Construct from components
        blastP1(const volScalarField& T);

        //- Construct from components
        blastP1(const dictionary& dict, const volScalarField& T);

        //- Disallow default bitwise copy construction
        blastP1(const blastP1&) = delete;


    //- Destructor
    virtual ~blastP1();


    // Member Functions

        // Edit

            //- Solve radiation equation(s)
            virtual void calculate();

            //- Read radiation properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (for power of T^4)
            virtual scalar cellRp(const label celli) const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;

            //- Source term component (constant)
            virtual scalar cellRu(const label celli) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blastP1&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
