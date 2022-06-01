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
    Foam::reactingSolidBlastThermo

Description
    Reacting, single phase solid thermo model

SourceFiles
    reactingSolidBlastThermo.C
    reactingSolidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef reactingSolidBlastThermo_H
#define reactingSolidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentSolidBlastThermo.H"
#include "basicBlastChemistryModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class reactingSolidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class reactingSolidBlastThermo
:
    public multicomponentSolidBlastThermo<Thermo>
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
    reactingSolidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    reactingSolidBlastThermo
    (
        const HashPtrTable<typename Thermo::thermoType, word, string::hash>& thermoData,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~reactingSolidBlastThermo();


    // Member Functions

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
    #include "reactingSolidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
