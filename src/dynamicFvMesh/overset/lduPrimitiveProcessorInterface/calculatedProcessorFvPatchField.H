/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::calculatedProcessorFvPatchField

Group
    grpGenericBoundaryConditions

Description
    processorFvPatchField type bypassing fvPatch

    Used to temporarily add updateInterfaceMatrix capabilities to a matrix
    during overset solving. Supplies:
    - patchNeighbourField functionality (cached in *this as per
      processorFvPatchField)
    - initEvaluate/evaluate: caching of patchNeighbourField (see above)
    - initInterfaceMatrixUpdate etc: adding of neighbouring data

SourceFiles
    calculatedProcessorFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef calculatedProcessorFvPatchField_H
#define calculatedProcessorFvPatchField_H

#include "lduPrimitiveProcessorInterface.H"
#include "coupledFvPatchField.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class calculatedProcessorFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class calculatedProcessorFvPatchField
:
    public processorLduInterfaceField,
    public coupledFvPatchField<Type>
{
protected:

    // Private data

         //- Local reference cast into the interface
         const lduPrimitiveProcessorInterface& procInterface_;

        // Sending and receiving

            //- Send buffer
            mutable Field<Type> sendBuf_;

            //- Receive buffer
            mutable Field<Type> receiveBuf_;

            //- Scalar send buffer
            mutable scalarField scalarSendBuf_;

            //- Scalar receive buffer
            mutable scalarField scalarReceiveBuf_;

            //- Outstanding request
            mutable label outstandingSendRequest_;

            //- Outstanding request
            mutable label outstandingRecvRequest_;



    // Private Member Functions

        void addToInternalField
        (
            scalarField& result,
            const bool add,
            const scalarField& coeffs,
            const scalarField& vals
        ) const;


public:

    //- Runtime type information
    TypeName("calculatedProcessor");


    // Constructors

        //- Construct from patch and internal field
        calculatedProcessorFvPatchField
        (
            const lduInterface& interface,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct as copy
        calculatedProcessorFvPatchField
        (
            const calculatedProcessorFvPatchField<Type>&
        ) = delete;

        //- Construct as copy setting internal field reference
        calculatedProcessorFvPatchField
        (
            const calculatedProcessorFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new calculatedProcessorFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~calculatedProcessorFvPatchField() = default;


    // Member functions

        // processorLduInterfaceField implementation

            //- Return communicator used for comms
            virtual label comm() const
            {
                return procInterface_.comm();
            }

            //- Return processor number
            virtual int myProcNo() const
            {
                return procInterface_.myProcNo();
            }


            //- Return neighbour processor number
            virtual int neighbProcNo() const
            {
                return procInterface_.myProcNo();
            }

            //- Is the transform required
            virtual bool doTransform() const
            {
                return false;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return procInterface_.forwardT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }


        // Access

            //- Return true if this patch field is coupled.
            //  Our field supplies coefficients to the fvMatrix so
            //  should behave as a processorFvPatchField (in
            //  addBoundarySource it should not add to the source)
            virtual bool coupled() const
            {
                return Pstream::parRun();
            }

            //- Return neighbour field of internal field
            virtual tmp<Field<Type>> patchNeighbourField() const;


        // Evaluation functions

            //- Is all data available
            virtual bool ready() const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
//                 const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                Field<Type>& result,
//                 const bool add,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
            }

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& result,
                const bool add,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
            }

        //- Return null transform between processor patches
        virtual const transformer& transform() const
        {
            return transformer::null;
        }

        //- Write
        //virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "calculatedProcessorFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
