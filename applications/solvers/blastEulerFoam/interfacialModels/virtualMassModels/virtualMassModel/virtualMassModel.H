/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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
    Foam::virtualMassModel

Description

SourceFiles
    virtualMassModel.C
    virtualMassModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef virtualMassModel_H
#define virtualMassModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "BlendedInterfacialModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                         Class virtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class virtualMassModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("virtualMassModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            virtualMassModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair,
                const bool registerObject
            ),
            (dict, pair, registerObject)
        );


    // Static data members

        //- Force coefficient dimensions
        static const dimensionSet dimK;

    // Constructors

        //- Construct from a dictionary and a phase pair
        virtualMassModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~virtualMassModel();


    // Selectors

        static autoPtr<virtualMassModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Return the virtual mass coefficient
        virtual tmp<volScalarField> K() const = 0;

        // Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                      Class blendedVirtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class blendedVirtualMassModel
:
    public BlendedInterfacialModel<virtualMassModel>
{
public:

    typedef virtualMassModel modelType;

    // Constructors

        //- Inherit base class constructors
        using BlendedInterfacialModel<virtualMassModel>::BlendedInterfacialModel;


    // Member Functions

        //- Return the drag coefficient K
        tmp<volScalarField> K() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
