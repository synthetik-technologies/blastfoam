/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::NusseltNumberModels::duct

Description
    Nusselt number correlation of Dittus–Boelter

    \verbatim
        "Heat transfer in automobile radiators of the tubular type"
        Dittus, F.W., Boelter, L.M.K.,
        University of California press, Berkeley, Calif,
        1930
    \endverbatim

SourceFiles
    ductNusseltNumber.C

\*---------------------------------------------------------------------------*/

#ifndef ductNusseltNumber_H
#define ductNusseltNumber_H

#include "NusseltNumberModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace NusseltNumberModels
{

/*---------------------------------------------------------------------------*\
                          Class duct Declaration
\*---------------------------------------------------------------------------*/

class duct
:
    public NusseltNumberModel
{
public:

    //- Runtime type information
    TypeName("duct");


    // Constructors

        //- Construct from components
        duct
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~duct();


    // Member Functions

        //- Nusselt number
        virtual tmp<volScalarField> Nu() const;

        //- Nusselt number for cell
        virtual scalar cellNu(const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace NusseltNumberModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
