/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::errorEstimators::fieldValue

Description
    Error is the value of the given field

SourceFiles
    fieldValue.C

\*---------------------------------------------------------------------------*/

#ifndef fieldValue_H
#define fieldValue_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "errorEstimator.H"
#include "NamedEnum.H"

namespace Foam
{
namespace errorEstimators
{

/*---------------------------------------------------------------------------*\
                           Class cellSize Declaration
\*---------------------------------------------------------------------------*/

class cellSize
:
    public errorEstimator
{
public:

    enum SizeType
    {
        VOLUME,
        CMPT,
        CHARACTERISTIC,
        MAG
    };
    static const NamedEnum<SizeType, 4> sizeTypeNames;


private:
// Private data

    //- Parameters
    SizeType sizeType_;
    labelList cmpts_;
    vector minDX_;
    vector maxDX_;

    //- Read components
    labelList readCmpts(Istream& is) const;

public:

    //- Runtime type information
    TypeName("cellSize");

    // Constructor
    cellSize(const fvMesh& mesh, const dictionary& dict, const word& name);


    //- Destructor
    virtual ~cellSize();


    // Member Functions

        //- Update error
        virtual void update(const bool scale = true);

        //- Return non constant reference to error field
        virtual labelList maxRefinement() const;

        //- Read parameters
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorEstimators
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
