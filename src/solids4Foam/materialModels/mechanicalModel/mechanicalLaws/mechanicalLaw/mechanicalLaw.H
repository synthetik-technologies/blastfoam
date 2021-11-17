/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mechanicalLaw

Description
    Mechanical constitutive law for the solidModels.

SourceFiles
    mechanicalLaw.C
    newMechanicalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalLaw_H
#define mechanicalLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mechanicalLaw Declaration
\*---------------------------------------------------------------------------*/

class mechanicalLaw
{
    // Private data

        //- Name
        const word name_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Mechanical law dictionary
        const dictionary dict_;

        //- Name of the baseMesh region
        word baseMeshRegionName_;

        //- Nonlinear geometry enumerator
        const nonLinearGeometry::nonLinearType nonLinGeom_;

        //- Reference to density field
        const volScalarField& rho_;

        // Total deformation gradient volField
        mutable autoPtr<volTensorField> FPtr_;

        // Total deformation gradient surfaceField
        mutable autoPtr<surfaceTensorField> FfPtr_;

        //- Relative deformation gradient volField
        mutable autoPtr<volTensorField> relFPtr_;

        //- Relative deformation gradient surfaceField
        mutable autoPtr<surfaceTensorField> relFfPtr_;

        //- Jacobian of the toal deformation gradient volField
        mutable autoPtr<volScalarField> JPtr_;

        //- Jacobian of the toal deformation gradient surfaceField
        mutable autoPtr<surfaceScalarField> JfPtr_;

        // Hydrostatic stress (negative of hydrostatic pressure) volField
        autoPtr<volScalarField> sigmaHydPtr_;

        // Gradient of the hydrostatic stress field
        autoPtr<volVectorField> gradSigmaHydPtr_;

        //- Is the deformation of the solid model used
        bool useSolidDeformation_;


    // Private Member Functions

        //- Make the F field
        void makeF() const;

        //- Make the Ff field
        void makeFf() const;

        //- Make the relF field
        void makeRelF() const;

        //- Make the relFf field
        void makeRelFf() const;

        //- Make the J field
        void makeJ() const;

        //- Make the Jf field
        void makeJf() const;

        //- Disallow copy construct
        mechanicalLaw(const mechanicalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalLaw&);


protected:

        //- Switch to enable solution of pressure Laplacian equation
        //  This can help quell oscillations in the hydrostatic stress
        const Switch solvePressureEqn_;

        //- Optional: it is possible to scale the amount of smoothing in the
        //  pressure equation with this coefficient
        const scalar pressureSmoothingScaleFactor_;

        //- Return const reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return const reference to dict
        const dictionary& dict() const
        {
            return dict_;
        }

        //- If the case is plane stress
        bool planeStress() const;

        //- Return the base mesh region name
        const word& baseMeshRegionName() const
        {
            return baseMeshRegionName_;
        }

        //- Return const access to the base mesh
        const fvMesh& baseMesh() const
        {
            return mesh().time().lookupObject<fvMesh>
            (
                baseMeshRegionName()
            );
        }

//         //- Read a field from the starting directory, or construct
//         template<class Type, class Mesh, class Patch>
//         tmp<GeometricField<Type, Mesh, Patch>> readIfPresent
//         (
//             const fvMesh& mesh,
//             const word& name,
//             const dimensioned<Type>& initVal
//         ) const;


        //- Return a reference to the F field
        const volTensorField& F() const;

        //- Return a reference to the Ff field
        const surfaceTensorField& Ff() const;

        //- Return a reference to the relF field
        const volTensorField& relF() const;

        //- Return a reference to the relF field
        volTensorField& relF();

        //- Return a reference to the relFf field
        const surfaceTensorField& relFf() const;

        //- Return a reference to the relFf field
        surfaceTensorField& relFf();

        //- Return a reference to the J field
        const volScalarField& J() const;

        //- Return a reference to the Jf field
        const surfaceScalarField& Jf() const;

        //- Return a reference to the sigmaHyd field
        volScalarField& sigmaHyd();

        //- Return a reference to the gradSigmaHyd field
        volVectorField& gradSigmaHyd();

        //- Update the deformation gradient
        //  The lienarised shear and bulk modulii are used to enforce linear
        //  elasticity in the case that the enforceLinear flag is tripped by the
        //  solid model
        //  true is returned if enforceLinear is true (linear elasticity is
        //  enforced)
        bool updateF
        (
            volSymmTensorField& sigma,
            const dimensionedScalar& mu,
            const dimensionedScalar& K
        );

        //- Equivalent to the updateF function, except instead for the Ff
        //  surface field
        bool updateFf
        (
            surfaceSymmTensorField& sigma,
            const dimensionedScalar& mu,
            const dimensionedScalar& K
        );

        //- Update the hydrostatic stress field by solving a diffusion equation.
        //  An explicit calculation of sigmaHyd is passed as an argument.
        //  Passing a zero will enforce incompressibility, although for
        //  covergence it may be better to pass a penalty calculation of
        //  sigmaHyd
        void updateSigmaHyd
        (
            const volScalarField& sigmaHydExplicit,
            const dimensionedScalar& impK
        );

        //- Lookup the enforceLinear Switch in the solidModel
        const Switch& enforceLinear() const;

        //- Does the solver use an incremental approach
        //  i.e. does it solve for DD as opposed to D
        bool incremental() const;

        //- Return nonlinear geometry enumerator
        nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinGeom_;
        }

public:

    //- Runtime type information
    TypeName("mechanicalLaw");


    // Declare run-time constructor selection table

        //- Mechanical law for linear geometry i.e. small strains and rotations
        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalLaw,
            linGeomMechLaw,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict,
                const nonLinearGeometry::nonLinearType& nonLinGeom
            ),
            (name, mesh, dict, nonLinGeom)
        );

        //- Mechanical law for nonlinear geometry i.e. finite strains and
        //  rotations
        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalLaw,
            nonLinGeomMechLaw,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict,
                const nonLinearGeometry::nonLinearType& nonLinGeom
            ),
            (name, mesh, dict, nonLinGeom)
        );


    // Selectors

        //- Create a linear geometry mechanical law
        static autoPtr<mechanicalLaw> NewLinGeomMechLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );

        //- Create a nonlinear geometry mechanical law
        static autoPtr<mechanicalLaw> NewNonLinGeomMechLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Constructors

        //- Construct from dictionary
        mechanicalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~mechanicalLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return density
        virtual tmp<volScalarField> rho() const
        {
            return rho_;
        }

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const = 0;

        //- Return the implicit stiffness surface field
        //  This is the diffusivity for the Laplacian term
        virtual tmp<surfaceScalarField> impKf() const;

        //- Return the implicit stiffness for a patch
        //  This is the diffusivity for the Laplacian term
        virtual tmp<scalarField> impK(const label) const = 0;

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const = 0;

        //- Return elastic modulus
        virtual tmp<volScalarField> elasticModulus() const = 0;

        //- Return shear modulus
        virtual tmp<volScalarField> shearModulus() const = 0;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma) = 0;

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update total accumulated fields
        virtual void updateTotalFields()
        {}

        //- Return the desired new time-step
        virtual scalar newDeltaT();

        void setUseSolidDeformation()
        {
            useSolidDeformation_ = true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
