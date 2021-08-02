/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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
    Foam::solidBlastThermo

Description
    Basic solid class for returning thermal conductivity (isotropic and
    anisotropic)

SourceFiles
    solidBlastThermo.C
    newSolidBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidBlastThermo_H
#define solidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicBlastThermo.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class basicFluidThermo Declaration
\*---------------------------------------------------------------------------*/

class solidBlastThermo
:
    public basicBlastThermo
{

public:

    //- Runtime type information
    TypeName("solidBlastThermo");

    declareRunTimeSelectionTable
    (
        autoPtr,
        solidBlastThermo,
        basicSolid,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master,
            const word& masterName
        ),
        (name, mesh, dict, master, masterName)
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        solidBlastThermo,
        detonatingSolid,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master,
            const word& masterName
        ),
        (name, mesh, dict, master, masterName)
    );

    // Constructor

        //- Construct with input field
        solidBlastThermo
        (
            const word& phaseName,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        //- Construct and add fields to db
        solidBlastThermo
        (
            const word& phaseName,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

    //- Selector
        static autoPtr<solidBlastThermo> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        static autoPtr<solidBlastThermo> NewBasic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        static autoPtr<solidBlastThermo> NewDetonating
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

    //- Destructor
    virtual ~solidBlastThermo();


        //- Return energy source
        virtual tmp<volScalarField> ESource() const = 0;

        //- Thermal conductivity [W/m/K]
        virtual tmp<volVectorField> Kappa() const = 0;

        //- Anisotropic thermal conductivity [W/m/K]
        virtual tmp<vectorField> Kappa
        (
            const label patchi
        ) const = 0;

        //- Dynamic viscosity of mixture [kg/m/s]
        virtual tmp<volScalarField> mu() const;

        //- Dynamic viscosity of mixture for patch [kg/m/s]
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Kinematic viscosity of mixture [m^2/s]
        virtual tmp<volScalarField> nu() const;

        //- Kinematic viscosity of mixture for patch [m^2/s]
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Kinematic viscosity of mixture for cell [m^2/s]
        virtual scalar nui(const label celli) const;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
