/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    mechanicalEnergies

Description
    Class used to calculate the energies involved in a solid mechanics
    simulation, including the external work, internal strain energy, and
    kinetic energy.

    Energies are integrated using the trapezoidal rule.

    The class also calculates a linear bulk viscous pressure term, which is used
    to dissipate high frequency energy in explicit simulations.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    mechanicalEnergies.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalEnergies_H
#define mechanicalEnergies_H

#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class mechanicalEnergies Declaration
\*---------------------------------------------------------------------------*/

class mechanicalEnergies
{
    // Private data

        // Reference to mesh
        const fvMesh& mesh_;

        //- Accumulated external work
        scalar externalWork_;

        //- Accumulated external work old time
        scalar externalWorkOldTime_;

        //- Accumulated internal energy
        scalar internalEnergy_;

        //- Accumulated internal energy old time
        scalar internalEnergyOldTime_;

        //- Kinetic energy
        scalar kineticEnergy_;

        //- Accumulated energy dissipated due to Laplacian smoothing term
        // scalar laplacianSmoothingEnergy_;

        //- Accumulated energy dissipated due to Laplacian smoothing term old
        //  time
        // scalar laplacianSmoothingEnergyOldTime_;

        //- Linear bulk viscosity energy
        scalar linearBulkViscosityEnergy_;

        //- Linear bulk viscosity energy old time
        scalar linearBulkViscosityEnergyOldTime_;

        //- Linear bulk viscosity coefficient
        //  Defaults to 0.06, as in Abaqus
        const scalar linearBulkViscosityCoeff_;

        //- Bulk viscosity stablisation pressure
        autoPtr<surfaceScalarField> viscousPressurePtr_;

        //- Volumetric strain used calculate the linear bulk viscosity pressure
        autoPtr<volScalarField> epsilonVolPtr_;

        //- Output file for the energies
        autoPtr<OFstream> energiesFilePtr_;

        //- Time index to know when a new time step occurs
        label curTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mechanicalEnergies(const mechanicalEnergies&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalEnergies&);

public:

    //- Runtime type information
    TypeName("mechanicalEnergies");

    // Constructors

        //- Construct from components
        mechanicalEnergies
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Destructor

        virtual ~mechanicalEnergies()
        {}


    // Member Functions

        //- Viscous pressure
        const surfaceScalarField& viscousPressure
        (
            const volScalarField& rho,
            const surfaceScalarField& waveSpeed,
            const volTensorField& gradD
        );

        //- Volumetric strain
        const volScalarField& epsilonVol
        (
            const volTensorField& gradD
        );

        //- Check global energy conservation
        void checkEnergies
        (
            const volScalarField& rho,
            const volVectorField& U,
            const volVectorField& D,
            const volVectorField& DD,
            const volSymmTensorField& sigma,
            const volTensorField& gradD,
            const volTensorField& gradDD,
            const surfaceScalarField& waveSpeed,
            const dimensionedVector& g,
            const scalar, // laplacianSmoothCoeff,
            const surfaceScalarField& impKf
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
