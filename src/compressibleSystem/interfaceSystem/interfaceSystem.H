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
    Foam::interfaceSystem

Description
    Class for general interface functions

SourceFiles
    interfaceSystem.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceSystem_H
#define interfaceSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvc.H"
#include "reconstruction.H"
#include "HashPtrTable.H"

namespace Foam
{
class surfaceTensionModel;

/*---------------------------------------------------------------------------*\
                Class interfaceSystem Declaration
\*---------------------------------------------------------------------------*/

class interfaceSystem
{
public:

    static wordHashSet compressionSchemes;

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

    typedef HashPtrTable<surfaceTensionModel, interfacePair, Hash<interfacePair>>
        surfaceTensionTable;


protected:
// Protected data

    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Velocity
    const volVectorField& velocity_;

    //- Stabilisation for normalisation of the interface normal
    dimensionedScalar deltaN_;

    //- Interface thicknesses
    HashTable<scalar> dxs_;


    // Protected Member Functions

        //- Return normalized interface normal vector
        tmp<volVectorField> nHat(const volScalarField& alpha) const;

        //- Return normalized interface normal vector
        tmp<volVectorField> nHat
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Return normalized interpolated interface normal vector
        tmp<surfaceVectorField> nHatfv(const volScalarField& alpha) const;

        //- Return normalized interpolated interface normal vector
        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Return interface normal flux
        tmp<surfaceScalarField> nHatf(const volScalarField& alpha) const;

        //- Return interface normal flux
        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Correct boundary interface normal
        void correctContactAngle
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            surfaceVectorField::Boundary& nHatb
        ) const;

        //- Return Curvature
        tmp<volScalarField> K
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Correct interpolated values so that upwinding at the interface uses
        //  The interface normal direction rather than the velocity
        //  This is only for phase quantities
        template<class Type>
        static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correctInterfaceField
        (
            const word& schemeKey,
            const GeometricField<Type, fvPatchField, volMesh>& alpha,
            const GeometricField<Type, fvPatchField, volMesh>& psi,
            const GeometricField<scalar, fvsPatchField, surfaceMesh>& nHatf,
            GeometricField<Type, fvsPatchField, surfaceMesh>& psif
        );

        template<class Type>
        static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correctInterfaceField
        (
            const word& schemeKey,
            const GeometricField<Type, fvPatchField, volMesh>& alpha,
            const GeometricField<Type, fvPatchField, volMesh>& psi,
            const tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>>& nHatf,
            GeometricField<Type, fvsPatchField, surfaceMesh>& psif
        );

        template<class Type>
        static void correctInterfaceField
        (
            const GeometricField<Type, fvPatchField, volMesh>& alpha,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& psiIntf,
            GeometricField<Type, fvsPatchField, surfaceMesh>& psif
        );

        template<class Type>
        static void correctInterfaceField
        (
            const GeometricField<Type, fvPatchField, volMesh>& alpha,
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& psiIntf,
            GeometricField<Type, fvsPatchField, surfaceMesh>& psif
        );

        template<class Type>
        static void correctInterfaceField
        (
            const word& schemeKey,
            const GeometricField<Type, fvPatchField, volMesh>& alpha,
            const GeometricField<Type, fvPatchField, volMesh>& psi,
            const GeometricField<scalar, fvsPatchField, surfaceMesh>& nHatf,
            GeometricField<Type, fvsPatchField, surfaceMesh>& psiOwn,
            GeometricField<Type, fvsPatchField, surfaceMesh>& psiNei
        );

        // Return interface thickness
        dimensionedScalar dx(const IOobject& io) const
        {
            if (dxs_.found(io.group()))
            {
                return dimensionedScalar("dx", dimLength, dxs_[io.group()]);
            }
            else
            {
                return dimensionedScalar("dx", dimLength, dxs_["default"]);
            }
        }

public:

    TypeName("multiphaseInterface");

    // Constructor
    interfaceSystem(const volVectorField& U, const dictionary& dict);


    //- Destructor
    virtual ~interfaceSystem();


    // Member functions

        bool usesCompression(const word& name) const;

        //- Update fluxes
        void updateDeltaN();

        tmp<surfaceScalarField> surfaceTensionForce
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const surfaceTensionModel& stModel
        ) const;

        tmp<surfaceScalarField> surfaceTensionForce
        (
            const UPtrList<volScalarField>& alphas,
            const surfaceTensionTable& stModels
        ) const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface
        (
            const volScalarField& alpha
        ) const;

        tmp<volScalarField> nearInterface
        (
            const volScalarField& psi,
            const dimensionedScalar& dx
        ) const;

        tmp<surfaceScalarField> nearInterface
        (
            const surfaceScalarField& alphaf
        ) const;

        //- Return a level-set function determined from a volume fraction
        tmp<volScalarField> levelSet
        (
            const volScalarField& alpha
        ) const;

        //- redistance the level set function
        void redistance(volScalarField& psi) const;

        //- Using the ghost fluid method, extend the given field
        void GFMExtension
        (
            volScalarField& x,
            const volScalarField& psi
        ) const;

        void GFMExtension
        (
            volScalarField& x,
            const volScalarField& psi,
            const surfaceScalarField& nHatf,
            const volScalarField& divnHatf
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interfaceSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
