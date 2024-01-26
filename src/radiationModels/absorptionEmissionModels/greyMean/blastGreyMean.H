/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::radiationModels::absorptionEmission::blastGreyMean

Description
    blastGreyMean radiation absorption and emission coefficients for continuous phase

    The coefficients for the species in the Look up table have to be specified
    for use in moles x P [atm], i.e. (k[i] = species[i]*p*9.869231e-6).

    The coefficients for CO and soot or any other added are multiplied by the
    respective mass fraction being solved

    All the species in the dictionary need either to be in the look-up table or
    being solved. Conversely, all the species solved do not need to be included
    in the calculation of the absorption coefficient

    The names of the species in the absorption dictionary must match exactly the
    name in the look-up table or the name of the field being solved

    The look-up table ("speciesTable") file should be in constant

Example
    \verbatim
    absorptionEmissionModel blastGreyMean;

    blastGreyMeanCoeffs
    {
        lookUpTableFileName     "speciesTable";

        EhrrCoeff       0.0;

        CO2
        {
            Tcommon     300.;   // Common Temp
            invTemp     true;   // Is the polynomial using inverse temperature?
            Tlow        300.;   // Low Temp
            Thigh       2500.;  // High Temp

            loTcoeffs           // coeffs for T < Tcommon
            (
                0               //  a0            +
                0               //  a1*T          +
                0               //  a2*T^(+/-)2   +
                0               //  a3*T^(+/-)3   +
                0               //  a4*T^(+/-)4   +
                0               //  a5*T^(+/-)5   +
            );
            hiTcoeffs           // coeffs for T > Tcommon
            (
                18.741
                -121.31e3
                273.5e6
                -194.05e9
                56.31e12
                -5.8169e15
            );
        }
    }
    \endverbatim

SourceFiles
    blastGreyMean.C

\*---------------------------------------------------------------------------*/

#ifndef blastGreyMean_H
#define blastGreyMean_H

#include "interpolationLookUpTable.H"
#include "blastAbsorptionEmissionModel.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "fluidThermo.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicSpecieMixture;

namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                          Class blastGreyMean Declaration
\*---------------------------------------------------------------------------*/

class blastGreyMean
:
    public blastAbsorptionEmissionModel
{
public:

    // Public data

        // Absorption Coefficients
        List<absorptionCoeffs> coeffs_;


protected:

    // Protected data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Hash table of species names
        HashTable<label> speciesNames_;

        //- Indices of species in the look-up table
        List<label> specieIndex_;

        //- Phase indices of species
        List<labelList> speciePhases_;

        //- Mol fractions (constant values)
        List<scalar> Xj_;

        //- Look-up table of species related to ft
        mutable autoPtr<interpolationLookUpTable> lookUpTablePtr_;

        //- Thermo package
        const fluidThermo& thermo_;

        //- Pointer list of species in the registry involved in the absorption
        UPtrList<volScalarField> Yj_;

        //- Pointer list of phases in the registry involved in the absorption
        UPtrList<volScalarField> volumeFractions_;

        //- Pointer to specie mixture
        UautoPtr<const basicSpecieMixture> mixture_;

        //- Pointer to ft field
        UautoPtr<const volScalarField> ft_;


public:

    //- Runtime type information
    TypeName("greyMean");


    // Constructors

        //- Construct from components
        blastGreyMean
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelName=typeName
        );


    //- Destructor
    virtual ~blastGreyMean();


    // Member Functions

        //- Absorption coefficient for continuous phase
            virtual tmp<volScalarField> aCont(const label bandI = 0) const;

            //- Absorption coefficient for continuous phase
            virtual scalar cellaCont
            (
                const label celli,
                const label bandI = 0
            ) const;

            //- Emission coefficient for continuous phase
            virtual tmp<volScalarField> eCont(const label bandI = 0) const;

            //- Emission coefficient for continuous phase
            virtual scalar celleCont
            (
                const label celli,
                const label bandI = 0
            ) const;

            //- Emission contribution for continuous phase
            virtual tmp<volScalarField> ECont(const label bandI = 0) const;

            //- Emission contribution for continuous phase
            virtual scalar cellECont
            (
                const label celli,
                const label bandI = 0
            ) const;

        //-
        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
