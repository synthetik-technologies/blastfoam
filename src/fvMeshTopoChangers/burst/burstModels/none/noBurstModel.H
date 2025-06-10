/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
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
    Foam::burstModels::none

Description
    No bursting

SourceFiles
    noBurstModel.C
    noBurstModel.C

\*---------------------------------------------------------------------------*/

#ifndef noBurstModel_H
#define noBurstModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "burstModel.H"


namespace Foam
{
namespace burstModels
{

/*---------------------------------------------------------------------------*\
                           Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public burstModel
{
public:

    //- Runtime type information
    TypeName("none");

    //- Constructor
    none
    (
        const dictionary& dict,
        const bool coupled
    );

    //- Clone
    virtual autoPtr<burstModel> clone() const
    {
        return autoPtr<burstModel>(new none());
    }

    //- Destructor
    virtual ~none();


    // Member Functions

        //- Find and set face indices to change, return if any added
        virtual bool facesToChange
        (
            const fvPatch& patch1,
            const fvPatch& patch2,
            labelList& masterFaces,
            labelList& slaveFaces
        ) const;

        //- Return if burst occurs
        virtual label update
        (
            const objectRegistry& obr,
            const labelList& ownCells,
            const labelList& neiCells,
            const scalarField& V
        ) const;

        //- Write
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace burstModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
