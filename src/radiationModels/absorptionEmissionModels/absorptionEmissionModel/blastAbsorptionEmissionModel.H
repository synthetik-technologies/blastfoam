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
    Foam::radiationModels::blastAbsorptionEmissionModel

Description
    Model to supply absorption and emission coefficients for radiation
    modelling

\*---------------------------------------------------------------------------*/

#ifndef blastAbsorptionEmissionModel_H
#define blastAbsorptionEmissionModel_H

#include "absorptionEmissionModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                   Class blastAbsorptionEmissionModel Declaration
\*---------------------------------------------------------------------------*/

class blastAbsorptionEmissionModel
:
    public absorptionEmissionModel
{
public:

    //- Runtime type information
    TypeName("blastAbsorptionEmissionModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            blastAbsorptionEmissionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        blastAbsorptionEmissionModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Selector
    static autoPtr<blastAbsorptionEmissionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~blastAbsorptionEmissionModel();


    // Member Functions

        // Access

            //- Reference to the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }


            // Absorption coefficient

                //- Absorption coefficient (net) for celli
                virtual scalar cella
                (
                    const label celli,
                    const label bandI = 0
                ) const;

                //- Absorption coefficient for continuous phase
                virtual scalar cellaCont
                (
                    const label celli,
                    const label bandI = 0
                ) const;

                //- Absorption coefficient for dispersed phase
                virtual scalar cellaDisp
                (
                    const label celli,
                    const label bandI = 0
                ) const;


            // Emission coefficient

                //- Emission coefficient (net)
                virtual scalar celle
                (
                    const label celli,
                    const label bandI = 0
                ) const;

                //- Return emission coefficient for continuous phase
                virtual scalar celleCont
                (
                    const label celli,
                    const label bandI = 0
                ) const;

                //- Return emission coefficient for dispersed phase
                virtual scalar celleDisp
                (
                    const label celli,
                    const label bandI = 0
                ) const;


            // Emission contribution

                //- Emission contribution (net)
                virtual scalar cellE
                (
                    const label celli,
                    const label bandI = 0
                ) const;

                //- Emission contribution for continuous phase
                virtual scalar cellECont
                (
                    const label celli,
                    const label bandI = 0
                ) const;

                //- Emission contribution for dispersed phase
                virtual scalar cellEDisp
                (
                    const label celli,
                    const label bandI = 0
                ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
