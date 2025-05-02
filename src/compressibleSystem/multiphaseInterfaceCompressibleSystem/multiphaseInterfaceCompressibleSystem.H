/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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
    Foam::multiphaseInterfaceCompressibleSystem

Description
    Compressible system containing more than two equation of states.

SourceFiles
    multiphaseInterfaceCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseInterfaceCompressibleSystem_H
#define multiphaseInterfaceCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiphaseCompressibleSystem.H"
#include "multiphaseFluidBlastThermo.H"
#include "surfaceTensionModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multiphaseInterfaceCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseInterfaceCompressibleSystem
:
    public multiphaseCompressibleSystem
{
public:

    class interfacePair
    :
        public Pair<word>
    {
    public:

        class hash
        :
            public Hash<interfacePair>
        {
        public:

            hash()
            {}

            label operator()(const interfacePair& key) const
            {
                return word::hash()(key.first()) + word::hash()(key.second());
            }
        };


        // Constructors

            using Pair<word>::Pair;

            interfacePair(const volScalarField& alpha1, const volScalarField& alpha2)
            :
                Pair<word>(alpha1.group(), alpha2.group())
            {}

            word name() const
            {
                return first() + "_" + second();
            }

        // Friend Operators

            friend bool operator==(const interfacePair& a, const interfacePair& b)
            {
                return
                (
                    ((a.first() == b.first()) && (a.second() == b.second()))
                 || ((a.first() == b.second()) && (a.second() == b.first()))
                );
            }

            friend bool operator!=(const interfacePair& a, const interfacePair& b)
            {
                return (!(a == b));
            }
    };

    typedef HashPtrTable<surfaceTensionModel, interfacePair, interfacePair::hash>
        surfaceTensionTable;

private:

    //- Surface tension models
    surfaceTensionTable stModels_;

    //- Stabilisation for normalisation of the interface normal
    dimensionedScalar deltaN_;

    // Protected member functions

        //- Return Curvature
        tmp<volScalarField> K
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Add source terms to momentum and energy equations
        virtual void addSources
        (
            volVectorField::Internal& rhoUSource,
            volScalarField::Internal& rhoESource
        ) const;


public:

    TypeName("multiphaseInterface");

    // Constructor
    multiphaseInterfaceCompressibleSystem
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const bool initialize = true
    );


    //- Destructor
    ~multiphaseInterfaceCompressibleSystem();


    // Member functions

        //- Update fluxes
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
