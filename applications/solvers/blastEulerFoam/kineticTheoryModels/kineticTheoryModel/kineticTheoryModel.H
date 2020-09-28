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
    Foam::kineticTheoryModel

Description
    Kinetic theory particle phase RAS model

    Reference:
    \verbatim
        van Wachem, B. G. M. (2000).
        Derivation, implementation, and validation of computer simulation models
        for gas-solid fluidized beds.
        PhD Thesis, TU Delft.
    \endverbatim

    There are no default model coefficients.

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "phaseModel.H"
#include "phasePairKey.H"
#include "dragModel.H"
#include "kineticTheorySystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
:
    public tmp<kineticTheoryModel>::refCount
{
private:
    // Private Data

            typedef
                HashTable
                <
                    autoPtr<phasePair>,
                    phasePairKey,
                    phasePairKey::hash
                >
                phasePairTable;

            //- Const refrence to phase
            const phaseModel& phase_;


protected:
    // Protected data

        // Kinetic Theory Model coefficients

            //- Pointer to central kineticTheory model
            kineticTheorySystem& kineticTheorySystem_;

            //- Maximum turbulent viscosity
            dimensionedScalar maxNut_;


        // Kinetic Theory Model Fields

            //- The granular energy/temperature
            volScalarField Theta_;

            //- Radial distribution function for single particle phase
            volScalarField gs0_;

            //- Radial distribution function for single particle phase
            volScalarField gs0Prime_;

            //- The granular bulk viscosity
            volScalarField lambda_;

            //- Granular pressure
            volScalarField Ps_;

            //- Frictional pressure
            volScalarField Pfric_;

            //- Total pressure
            volScalarField Ptot_;

            //- The granular "thermal" conductivity
            volScalarField kappa_;

            //- Viscosity
            volScalarField nut_;

            //- Frictional viscosity
            volScalarField nuFric_;

            //- Total viscosity
            volScalarField nuTotal_;

            //- Coefficient of restitution
            scalar es_;



    // Private Member Functions

        void correctNut()
        {}

        //- Lookup or construct kineticTheorySystem
        kineticTheorySystem& lookupOrConstruct(const fvMesh& mesh, const char*) const;


public:

    //- Runtime type information
    TypeName("kineticTheory");


    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const phaseModel& phase,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        kineticTheoryModel(const kineticTheoryModel&) = delete;


    //- Destructor
    virtual ~kineticTheoryModel();


    // Member Functions

        //- Return true if viscosity is included
        bool includeViscosity() const
        {
            return kineticTheorySystem_.includeViscosity();
        }

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const
        {
            return nut_.boundaryField()[patchi];
        }

        tmp<volScalarField> nuEff() const
        {
            return nuTotal_;
        }

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Solve the kinetic theory equations and correct the viscosity
        virtual void correct();

        //- Return cont reference to phase model
        const phaseModel& phase() const;

        //- Return the phase packing limit
        scalar alphaMax() const
        {
            return phase_.alphaMax();
        }

        //- Return const reference to granular temperature field
        const volScalarField& Theta() const;

        //- Return non-const reference to granular temperature field
        volScalarField& Theta();

        //- Return const reference to radial distribution function
        const volScalarField& gs0() const;

        //- Return const reference to radial distribution function
        const volScalarField& gs0Prime() const;

        //- Return coefficient of restitution
        const scalar& es() const
        {
            return es_;
        }

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
