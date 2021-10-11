/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::dynamicOversetBlastFvMesh

Description
    dynamicFvMesh with support for overset meshes.

SourceFiles
    dynamicOversetBlastFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicOversetBlastFvMesh_H
#define dynamicOversetBlastFvMesh_H

#include "dynamicMotionSolverListBlastFvMesh.H"
#include "labelIOList.H"
#include "fvMeshPrimitiveLduAddressing.H"
#include "lduInterfaceFieldPtrsList.H"
#include "volFields.H"
// #include "PackedBoolList.H"
#include "SolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistribute;
class lduPrimitiveProcessorInterface;
class GAMGAgglomeration;

/*---------------------------------------------------------------------------*\
                   Class dynamicOversetBlastFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicOversetBlastFvMesh
:
    public dynamicMotionSolverListBlastFvMesh
{
protected:

    // Protected data

        //- Select base addressing (false) or locally stored extended
        //  lduAddressing (true)
        mutable bool active_;

        //- Extended addressing (extended with local interpolation stencils)
        mutable autoPtr<fvMeshPrimitiveLduAddressing> lduPtr_;

        //- Added (processor)lduInterfaces for remote bits of stencil.
        //PtrList<const lduInterface> remoteStencilInterfaces_;
        mutable PtrList<const lduPrimitiveProcessorInterface>
            remoteStencilInterfaces_;

        //- Interfaces for above mesh. Contains both original and
        //  above added processorLduInterfaces
        mutable lduInterfacePtrsList allInterfaces_;

        //- Corresponding faces (in above lduPtr) to the stencil
        mutable labelListList stencilFaces_;

        //- Corresponding patches (in above lduPtr) to the stencil
        mutable labelListList stencilPatches_;

        //- From old to new face labels
        mutable labelList reverseFaceMap_;


   // Protected Member Functions

        //- Calculate the extended lduAddressing
        virtual bool updateAddressing() const;

        //- Debug: print matrix
        template<class Type>
        void write
        (
            Ostream&,
            const fvMatrix<Type>&,
            const lduAddressing&,
            const lduInterfacePtrsList&
        ) const;

        //- Explicit interpolation of acceptor cells from donor cells
        template<class T>
        void interpolate(Field<T>& psi) const;

        //- Explicit interpolation of acceptor cells from donor cells with
        //  boundary condition handling
        template<class GeoField>
        void interpolate(GeoField& psi) const;

        //- Helper: strip off trailing _0
        static word baseName(const word& name);

        //- Explicit interpolation of all registered fields. Excludes
        //  selected fields (and their old-time fields)
        template<class GeoField>
        void interpolate(const wordHashSet& suppressed);

        //- Freeze values at holes
        //template<class Type>
        //void freezeHoles(fvMatrix<Type>&) const;

        //- Get scalar interfaces of certain type
        //template<class GeoField, class PatchType>
        //lduInterfaceFieldPtrsList scalarInterfaces(const GeoField& psi) const;

        //- Determine normalisation for interpolation. This equals the
        //  original diagonal except stabilised for zero diagonals (possible
        //  in hole cells)
        template<class Type>
        tmp<scalarField> normalisation(const fvMatrix<Type>& m) const;

        //- Add interpolation to matrix (coefficients)
        template<class Type>
        void addInterpolation(fvMatrix<Type>&, const scalarField& norm) const;

        //- Solve given dictionary with settings
        template<class Type>
        SolverPerformance<Type> solve(fvMatrix<Type>&, const dictionary&) const;

        //- Debug: correct coupled bc
        template<class GeoField>
        static void correctCoupledBoundaryConditions(GeoField& fld);

        //- Average norm of valid neighbours
        scalar cellAverage
        (
            const labelList& types,
            const labelList& nbrTypes,
            const scalarField& norm,
            const scalarField& nbrNorm,
            const label celli,
            PackedBoolList& isFront
        ) const;

        //- Debug: dump agglomeration
        void writeAgglomeration
        (
            const GAMGAgglomeration& agglom
        ) const;


private:

    // Private Member Functions

        //- No copy construct
        dynamicOversetBlastFvMesh(const dynamicOversetBlastFvMesh&) = delete;

        //- No copy assignment
        void operator=(const dynamicOversetBlastFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicOversetBlastFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicOversetBlastFvMesh(const IOobject& io);//, const bool doInit=true);


    //- Destructor
    virtual ~dynamicOversetBlastFvMesh();


    // Member Functions


        // Extended addressing

            //- Return extended ldu addressing
            const fvMeshPrimitiveLduAddressing& primitiveLduAddr() const;

            //- Return ldu addressing. If active: is (extended)
            //  primitiveLduAddr
            virtual const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set. If active:
            //  return additional remoteStencilInterfaces_
            virtual lduInterfacePtrsList interfaces() const;

            //- Return old to new face addressing
            const labelList& reverseFaceMap() const
            {
                return reverseFaceMap_;
            }

            //- Return true if using extended addressing
            bool active() const
            {
                return active_;
            }

            //- Enable/disable extended addressing
            void active(const bool f) const
            {
                active_ = f;

                if (active_)
                {
                    DebugInfo<< "Switching to extended addressing with nFaces:"
                        << primitiveLduAddr().lowerAddr().size()
                        << endl;
                }
                else
                {
                    DebugInfo<< "Switching to base addressing with nFaces:"
                        << fvMesh::lduAddr().lowerAddr().size()
                        << endl;
                }
            }


        // Overset

            // Explicit interpolation

                virtual void interpolate(scalarField& psi) const
                {
                    interpolate<scalar>(psi);
                }

                virtual void interpolate(vectorField& psi) const
                {
                    interpolate<vector>(psi);
                }

                virtual void interpolate(sphericalTensorField& psi) const
                {
                    interpolate<sphericalTensor>(psi);
                }

                virtual void interpolate(symmTensorField& psi) const
                {
                    interpolate<symmTensor>(psi);
                }

                virtual void interpolate(tensorField& psi) const
                {
                    interpolate<tensor>(psi);
                }

                virtual void interpolate(volScalarField& psi) const
                {
                    interpolate<volScalarField>(psi);
                }

                virtual void interpolate(volVectorField& psi) const
                {
                    interpolate<volVectorField>(psi);
                }

                virtual void interpolate(volSphericalTensorField& psi) const
                {
                    interpolate<volSphericalTensorField>(psi);
                }

                virtual void interpolate(volSymmTensorField& psi) const
                {
                    interpolate<volSymmTensorField>(psi);
                }

                virtual void interpolate(volTensorField& psi) const
                {
                    interpolate<volTensorField>(psi);
                }


            // Implicit interpolation (matrix manipulation)

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<scalar> solve
                (
                    fvMatrix<scalar>& m,
                    const dictionary& dict
                ) const
                {
                    return solve<scalar>(m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<vector> solve
                (
                    fvMatrix<vector>& m,
                    const dictionary& dict
                ) const
                {
                    return solve<vector>(m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<symmTensor> solve
                (
                    fvMatrix<symmTensor>& m,
                    const dictionary& dict
                ) const
                {
                    return solve<symmTensor>(m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<tensor> solve
                (
                    fvMatrix<tensor>& m,
                    const dictionary& dict
                ) const
                {
                    return solve<tensor>(m, dict);
                }


        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update fields when mesh is updated
        virtual bool interpolateFields();

        //- Write using stream options
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write = true
        ) const;

        //- Debug: check halo swap is ok
        template<class GeoField>
        static void checkCoupledBC(const GeoField& fld);

        //- Correct boundary conditions of certain type (typeOnly = true)
        //  or explicitly not of the type (typeOnly = false)
        template<class GeoField, class PatchType>
        static void correctBoundaryConditions
        (
            typename GeoField::Boundary& bfld,
            const bool typeOnly
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "dynamicOversetBlastFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
