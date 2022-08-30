/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::EulerImplicitBlastChemistySolver

Description
    An Euler implicit solver for chemistry

    Euler implicit integration based on the reaction rate Jacobian is used to
    solver for the composition and Euler explicit integration for the
    temperature.  Separating the integration of temperature from composition
    in this manner is significantly more stable for exothermic systems

SourceFiles
    EulerImplicitBlastChemistySolver.C

\*---------------------------------------------------------------------------*/

#ifndef EulerImplicitBlastChemistySolver_H
#define EulerImplicitBlastChemistySolver_H

#include "blastChemistrySolver.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class EulerImplicitBlastChemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class EulerImplicitBlastChemistrySolver
:
    public blastChemistrySolver<ChemistryModel>
{
    // Private Data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Chemistry timescale coefficient
        scalar cTauChem_;

        //- Field encapsulating the composition, temperature and pressure
        mutable scalarField cTp_;

        //- Reaction rate field
        mutable scalarField R_;

        //- Reaction Jacobian
        mutable scalarSquareMatrix J_;

        //- Euler implicit integration matrix for composition
        mutable simpleMatrix<scalar> E_;


public:

    //- Runtime type information
    TypeName("EulerImplicit");


    // Constructors

        //- Construct from thermo
        EulerImplicitBlastChemistrySolver(const blastThermo& thermo);


    //- Destructor
    virtual ~EulerImplicitBlastChemistrySolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalar& p,
            scalar& T,
            scalarField& c,
            const label li,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EulerImplicitBlastChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
