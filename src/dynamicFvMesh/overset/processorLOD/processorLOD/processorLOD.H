/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::processorLOD

Description
    Base class to generate a parallel distribution map for sending sufficient
    target objects to cover a description of the source object, based on
    processor Level Of Detail (LOD) shapes

\*---------------------------------------------------------------------------*/

#ifndef processorLOD_H
#define processorLOD_H

#include "autoPtr.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistribute;

/*---------------------------------------------------------------------------*\
                        Class processorLOD Declaration
\*---------------------------------------------------------------------------*/

class processorLOD
{

protected:

    // Protected data

        //- Maximum number of objects per leaf
        label maxObjectsPerLeaf_;

        //- Number of objects of this type, e.g. number of faces/cells on this
        //- processor
        label nObjectsOfType_;


public:

    //- Runtime type information
    TypeName("processorLOD");

    //- Construct from components
    processorLOD
    (
        const label maxObjectsPerLeaf,
        const label nObjectsOfType
    );

    //- Destructor
    virtual ~processorLOD() = default;


    // Member Functions

        //- Return the parallel distribution map
        virtual autoPtr<mapDistribute> map() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
