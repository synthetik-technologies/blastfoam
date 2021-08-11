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
    Foam::phaseSolidBlastThermo

Description
    General fluid class that contains functions to return pressure and speed
    of sound

SourceFiles
    phaseSolidBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSolidBlastThermo_H
#define phaseSolidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicBlastThermo.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseSolidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class phaseSolidBlastThermo
:
    public basicBlastThermo
{
public:

    //- Runtime type information
    TypeName("solidBlastThermo");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseSolidBlastThermo,
            dictionary,
            (
                const word& phaseName,
                volScalarField& rho,
                volScalarField& e,
                volScalarField& T,
                const dictionary& dict,
                const word& masterName
            ),
            (phaseName, rho, e, T, dict, masterName)
        );


    // Constructor

        //- Construct with input field
        phaseSolidBlastThermo
        (
            const word& phaseName,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const word& masterName = word::null
        );


    // Selectors

        static autoPtr<phaseSolidBlastThermo> New
        (
            const word& name,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const word& masterName = word::null
        );

    //- Destructor
    virtual ~phaseSolidBlastThermo();


    //- Thermodynamic and transport functions

        // Access to transport state variables

            //- Calculate density
            virtual void updateRho() = 0;

            //- Calculate internal energy given a pressure
            virtual tmp<volScalarField> calce() const = 0;

            //- Return true if thermal conductivity is isotropic
            virtual bool isotropic() const = 0;

            //- Isotropic thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappa() const = 0;

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<volVectorField> Kappa() const = 0;

            //- Anisotropic thermal conductivity for patch [W/m/K]
            virtual tmp<vectorField> Kappa(const label patchi) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
