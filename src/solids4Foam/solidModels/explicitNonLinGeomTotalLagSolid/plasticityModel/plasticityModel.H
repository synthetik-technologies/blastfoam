/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::plasticityModel

Description
    Defines the plasticity model.

SourceFiles
    plasticityModel.C

\*---------------------------------------------------------------------------*/

#ifndef plasticityModel_H
#define plasticityModel_H

#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasticityModel
{
    // Private data

        //- Elastic left strain tensor
        volTensorField b_;

        //- First Piola Kirchhoff stress tensor
        volTensorField P_;

        //- First Piola Kirchhoff stress tensor
        volTensorField tau_;

        //- Von Mises stresses
        volScalarField vMises_;

        //- Plastic strain
        volScalarField strain_p_;

        //- Inverse of plastic Right Cauchy Green strain tensor
        volTensorField CpInv_;

        //- Pressure
        volScalarField p_;

        //- Energy of algorithm
        volScalarField energyAlgorithm_;

        //- Constitutive model
        word model_;

        //- Density
        dimensionedScalar rho_;

        //- Young's modulus of Elasticity
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Lame coefficient
        dimensionedScalar mu_;

        //- Lame coefficient
        dimensionedScalar lambda_;

        //- Bulk modulus
        dimensionedScalar kappa_;

        //- Hardening modulus
        dimensionedScalar Hm_;

        //- Initial yield stress
        dimensionedScalar Ys0_;

        //- Current yield stress
        volScalarField Ys_;

        //- Pressure wave speed
        dimensionedScalar pWaveSpeed_;

        //- Shear wave speed
        dimensionedScalar sWaveSpeed_;


private:

    //- Disallow default bitwise copy construct
        plasticityModel(const plasticityModel&);

    //- Disallow default bitwise assignment
        void operator=(const plasticityModel&);

public:

    //- Runtime type information
    TypeName("plasticityModel");


    // Constructors

        //- Construct from field and dictionary
        plasticityModel
        (
            const volTensorField&,
            const dictionary&
        );


    //- Destructor
    virtual ~plasticityModel();


    // Member functions

        // Access

            const volTensorField& piola() const
            {
                return P_;
            }

            const volScalarField& pressure() const
            {
                return p_;
            }

            const volScalarField& vonMisesStress() const
            {
                return vMises_;
            }

            const volScalarField& plasticStrain() const
            {
                return strain_p_;
            }

           const volScalarField& energyAlgorithm() const
            {
                return energyAlgorithm_;
            }

            const dimensionedScalar& density() const
            {
                return rho_;
            }

            const volScalarField& yieldStress() const
            {
                return Ys_;
            }

            const dimensionedScalar& pWaveSpeed() const
            {
                return pWaveSpeed_;
            }

            const dimensionedScalar& sWaveSpeed() const
            {
                return sWaveSpeed_;
            }

        // Edit

            //- Correct material model
            virtual void correct();

        // Output

            //- Print material properties
            void printMaterialProperties();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
