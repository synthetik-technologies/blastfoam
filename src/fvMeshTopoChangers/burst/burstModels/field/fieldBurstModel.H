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
    Foam::burstModels::field

Description
    Multiple fields used to determine if burst occurs

SourceFiles
    fieldBurstModel.C

\*---------------------------------------------------------------------------*/

#ifndef fieldBurstModel_H
#define fieldBurstModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "burstModel.H"


namespace Foam
{
namespace burstModels
{

/*---------------------------------------------------------------------------*\
                           Class field Declaration
\*---------------------------------------------------------------------------*/

class field
:
    public burstModel
{
protected:
// Protected variables

    //- Values that the patch will burst
    HashTable<scalar> burstValues_;


public:

    //- Runtime type information
    TypeName("field");

    //- Constructors

        //- Construct from dictionary
        field
        (
            const dictionary& dict,
            const bool coupled
        );

        //- Construct from dictionay and names/values
        field
        (
            const dictionary& dict,
            const bool coupled,
            const wordList& names,
            const scalarList& vals
        );

    //- Clone
    virtual autoPtr<burstModel> clone() const
    {
        return autoPtr<burstModel>(new field(*this));
    }


    //- Destructor
    virtual ~field();


    // Member Functions

        //- Find and set face indices to change, return if any added
        virtual bool facesToChange
        (
            const fvPatch& patch1,
            const fvPatch& patch2,
            labelList& masterFaces,
            labelList& slaveFaces
        ) const;

        //- Find and set face indices to change, return if any added
        virtual bool facesToChange
        (
            const fvPatch& patch,
            labelList& faces
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
