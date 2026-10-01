/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     | Website:  https://openfoam.org
   \\  /    A nd           | Copyright (C) 2019-2026
    \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::fv::pressurePhaseChange

Description
    Mass transfer from a source phase to a multi-component target phase
    triggered across a selected pressure interval.

SourceFiles
    pressurePhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef pressurePhaseChange_H
#define pressurePhaseChange_H

#include "fvModel.H"
#include "dimensionedScalar.H"
#include "HashTable.H"
#include "List.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations in namespace Foam
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class pressurePhaseChange Declaration
\*---------------------------------------------------------------------------*/

class pressurePhaseChange
:
    public fvModel
{
    // Private Data

        //- Struct representing a target phase with its proportion and species distribution
        struct targetPhase
        {
            word name;
            scalar fraction;
            HashTable<scalar> targetFractions;
        };

        //- Name of the source phase
        word sourcePhaseName_;

        //- List of target phases
        List<targetPhase> targetPhases_;

        //- Lower pressure threshold [Pa]
        scalar pLower_;

        //- Upper pressure threshold [Pa]
        scalar pUpper_;

        //- Smearing transition range [Pa]
        scalar pRange_;

        //- Relaxation time scale [s]
        scalar tau_;

        //- Are fractions specified as mole fractions
        bool isMoleFraction_;

        //- Integrated mass transfer rate [kg/s]
        scalar massRate_;

        //- Energy offset between reactant and product phases [J/kg]
        dimensionedScalar eOffset_;


    // Private Member Functions

        //- Return the smearing factor Gamma(p)
        tmp<volScalarField> gamma() const;

        //- Return the energy offset field [J/kg]
        tmp<volScalarField> energyOffset() const;


public:

    //- Runtime type information
    TypeName("pressurePhaseChange");


    // Constructors

        //- Construct from components
        pressurePhaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~pressurePhaseChange() = default;


    // Member Functions

        //- Return the list of fields for which the fvModel adds source term
        //  to the transport equation
        virtual wordList addSupFields() const;

        //- Return true if the fvModel adds a source term to the given
        //  field's transport equation
        virtual bool addsSupToField(const word& fieldName) const;


        // Mesh motion and topology changes

            //- Update for mesh motion
            virtual bool movePoints()
            {
                return false;
            }

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&)
            {}

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&)
            {}

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&)
            {}


        // Sources

            //- Pre-update mesh (called at start of integrate)
            virtual void preUpdateMesh();

            //- Main correction loop
            virtual void correct();

            //- Add explicit and implicit source terms to the continuity equation
            virtual void addSup
            (
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add source to momentum equation
            virtual void addSup
            (
                const volVectorField& field,
                fvMatrix<vector>& eqn
            ) const;

            //- Add explicit and implicit source terms to compressible continuity equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
