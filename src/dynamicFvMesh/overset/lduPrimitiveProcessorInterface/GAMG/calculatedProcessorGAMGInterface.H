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
    Foam::calculatedProcessorGAMGInterface

Description
    GAMG agglomerated processor interface.

SourceFiles
    calculatedProcessorGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef calculatedProcessorGAMGInterface_H
#define calculatedProcessorGAMGInterface_H

#include "GAMGInterface.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    template<class T>
    using LabelPairMap = HashTable<T, labelPair, Foam::Hash<labelPair>>;

    typedef HashSet<labelPair, Foam::Hash<labelPair>> labelPairHashSet;

    typedef HashTable<label, labelPair, Foam::Hash<labelPair>> labelPairLookup;


/*---------------------------------------------------------------------------*\
              Class calculatedProcessorGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class calculatedProcessorGAMGInterface
:
    public GAMGInterface,
    public processorLduInterface
{
    // Private data

        //- Communicator to use for parallel communication
        const label comm_;

        //- My processor rank in communicator
        label myProcNo_;

        //- Neighbouring processor rank in communicator
        label neighbProcNo_;

        //- Transformation tensor
        tensorField forwardT_;

        //- Message tag used for sending
        int tag_;


    // Private Member Functions

        //- No copy construct
        calculatedProcessorGAMGInterface
        (
            const calculatedProcessorGAMGInterface&
        ) = delete;

        //- No copy assignment
        void operator=(const calculatedProcessorGAMGInterface&) = delete;


public:

    //- Runtime type information
    TypeName("calculatedProcessor");

    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour restrict addressing
        calculatedProcessorGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Construct from components
        calculatedProcessorGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const labelUList& faceCells,
            const labelUList& faceRestrictAddresssing,
            const label coarseComm,
            const label myProcNo,
            const label neighbProcNo,
            const tensorField& forwardT,
            const int tag
        );

        //- Construct from Istream
        calculatedProcessorGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    //- Destructor
    virtual ~calculatedProcessorGAMGInterface() = default;


    // Member Functions

        // Interface transfer functions

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


        //- Processor interface functions

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

            virtual const transformer& transform() const
            {
                return transformer::I;
            }

        // I/O

            //- Write to stream
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
