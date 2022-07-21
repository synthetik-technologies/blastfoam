/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
31-03-2022 Synthetik Applied Technologies:  Added blast thermo models
-------------------------------------------------------------------------------
License
    This file is a derivative work of foam-extend.

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

Class
    thermalModel

Description
    Thermal material properties for solids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    thermalModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalModel_H
#define thermalModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "blastThermo.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidSubMeshes;

/*---------------------------------------------------------------------------*\
                         Class thermalModel Declaration
\*---------------------------------------------------------------------------*/

class thermalModel
{
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Optional solidSubMeshes
        UautoPtr<solidSubMeshes> subsetMeshes_;

        //- Thermal law
        autoPtr<blastThermo> thermoPtr_;

        //- Sub thermo models
        PtrList<blastThermo> thermoModels_;


    // Private Member Functions

        //- Disallow copy construct
        thermalModel(const thermalModel&);

        //- Disallow default bitwise assignment
        void operator=(const thermalModel&);


public:

    //- Runtime type information
    TypeName("thermalModel");

    // Constructors

        //- Construct from dictionary
        thermalModel(const fvMesh& mesh, const bool isSolid = true);


    // Destructor

        virtual ~thermalModel();


    // Member Functions

        //- Return const reference to the mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return thermal model
        const blastThermo& thermo() const
        {
            return thermoPtr_();
        }

        //- Access thermal model
        blastThermo& thermo()
        {
            return thermoPtr_();
        }

        //- Correct thermo models
        void correct();

        //- Access density
        volScalarField& rho();

        //- Return density
        const volScalarField& rho() const;

        //- Return specific heat capacity
        tmp<volScalarField> C() const;

        //- Return thermal conductivity
        tmp<volScalarField> k() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
