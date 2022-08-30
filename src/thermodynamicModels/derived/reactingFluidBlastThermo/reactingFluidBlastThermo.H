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
    Foam::reactingFluidBlastThermo

Description
    Non reacting, single phase fluid thermo model

SourceFiles
    reactingFluidBlastThermo.C
    reactingFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef reactingFluidBlastThermo_H
#define reactingFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentFluidBlastThermo.H"
#include "basicBlastChemistryModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reactingFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class reactingFluidBlastThermo
:
    public multicomponentFluidBlastThermo<Thermo>
{
// Private member data

    //- Chemistry model
    autoPtr<basicBlastChemistryModel> chemistry_;

    //- Solve chemistry at all ode steps
    bool odeChemistry_;


public:

    //- Runtime type information
    TypeName("reacting");

    // Constructor
    reactingFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    reactingFluidBlastThermo
    (
        const HashPtrTable<typename Thermo::thermoType, word, string::hash>& thermoData,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~reactingFluidBlastThermo();


    // Member Functions

        //- Pre update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "reactingFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
