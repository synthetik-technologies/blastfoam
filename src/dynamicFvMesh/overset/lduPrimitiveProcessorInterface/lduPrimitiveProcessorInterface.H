/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::lduPrimitiveProcessorInterface

Description
    Concrete implementation of processor interface. Used to temporarily store
    settings.

SourceFiles
    lduPrimitiveProcessorInterface.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveProcessorInterface_H
#define lduPrimitiveProcessorInterface_H

#include "GAMGInterface.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class lduPrimitiveProcessorInterface Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveProcessorInterface
:
    public lduInterface,
    public processorLduInterface
{
    // Private data

        const labelList faceCells_;

        //- My processor rank in communicator
        label myProcNo_;

        //- Neighbouring processor rank in communicator
        label neighbProcNo_;

        //- Transformation tensor
        tensorField forwardT_;

        //- Message tag used for sending
        int tag_;

        //- Communicator to use for parallel communication
        const label comm_;


    // Private Member Functions

        //- No copy construct
        lduPrimitiveProcessorInterface
        (
            const lduPrimitiveProcessorInterface&
        ) = delete;

        //- No copy assignment
        void operator=(const lduPrimitiveProcessorInterface&) = delete;


public:

    //- Runtime type information
    TypeName("calculatedProcessor");


    // Constructors

        //- Construct from components
        lduPrimitiveProcessorInterface
        (
            const labelUList& faceCells,
            const label myProcNo,
            const label neighbProcNo,
            const tensorField& forwardT,
            const int tag,
            const label comm = UPstream::worldComm
        );


    //- Destructor
    virtual ~lduPrimitiveProcessorInterface() = default;


    // Member Functions

        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return faceCells_;
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        // Processor interface functions

            //- Return communicator used for sending
            virtual label comm() const
            {
                return comm_;
            }

            //- Return processor number (rank in communicator)
            virtual int myProcNo() const
            {
                return myProcNo_;
            }

            //- Return neighbour processor number (rank in communicator)
            virtual int neighbProcNo() const
            {
                return neighbProcNo_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return forwardT_;
            }

            //- Return message tag used for sending
            virtual int tag() const
            {
                return tag_;
            }

            //- Return null transform between processor patches
            virtual const transformer& transform() const
            {
                return transformer::null;
            }

        // Edit

            //- Return message tag used for sending
            int& tag()
            {
                return tag_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
