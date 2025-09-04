/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
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
    Foam::diameterModels::none

Description
    Mass transfer rate based on arrhenius rate reaction for

    \verbatim
        Gurevich, M.A., Ozerova, G.E., Stepanov, A.M. (1970).
        Heterogeneous ignition of an aluminum particle in oxygen and water
        vapor.
        Combust Explos Shock Waves, 6 (3), 29-297.
    \endverbatim

SourceFiles
    noDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef noDiameterModel_H
#define noDiameterModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diameterModel.H"

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                    Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public diameterModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from runTime and dictionary
        none
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~none();



    // Public functions

        //- Return effective diameter
        virtual const volScalarField& d() const
        {
            NotImplemented;
            return NullObjectRef<volScalarField>();
        }

        //- Return surface area
        virtual tmp<volScalarField> A() const
        {
            NotImplemented;
            return tmp<volScalarField>();
        }

        //- Return volume
        virtual tmp<volScalarField> V() const
        {
            NotImplemented;
            return tmp<volScalarField>();
        }

        //- Return derivative of volume w.r.t diameter
        virtual tmp<volScalarField> dVdD() const
        {
            NotImplemented;
            return tmp<volScalarField>();
        }

        //- Return change in radius per unit time
        //  Only used in mass transfer
        virtual tmp<volScalarField> dMdt() const
        {
            NotImplemented;
            return tmp<volScalarField>();
        }

        //- Return the length to diameter ratio
        virtual tmp<volScalarField> LDratio() const
        {
            NotImplemented;
            return tmp<volScalarField>();
        }

        //- Return the length to diameter ratio for a cell
        virtual scalar LDratio(const label celli) const
        {
            NotImplemented;
            return 1.0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
