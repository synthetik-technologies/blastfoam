/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::boundaryTractions::hydrostaticPressure

Description
    Base class for constraining nodal displacement and velocity

SourceFiles
    displacementConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticPressureBoundaryTraction_H
#define hydrostaticPressureBoundaryTraction_H

#include "boundaryTraction.H"

namespace Foam
{
namespace boundaryTractions
{

class hydrostaticPressure
:
    public boundaryTraction
{
protected:

    //- Reference pressure
    const scalar pRef_;

    //- Reference height
    const scalar hRef_;

    //- Density of surrounding fluid
    const scalar rho_;

    //- Gravitational acceleration vector
    const vector g_;

    //- Magnitude of gravitational acceleration vector
    const scalar magg_;

    //- Copy constructor
    hydrostaticPressure(const hydrostaticPressure& stbt);

public:

    //- Runtime type information
    TypeName("hydrostaticPressure")

        // Constructor
        hydrostaticPressure
        (
            const word& name,
            const dictionary& dict,
            const feMesh1& femesh,
            const pointVectorField* DPtr
        );

        //- Clone
        virtual autoPtr<boundaryTraction> clone() const
        {
            return autoPtr<boundaryTraction>
            (
                new hydrostaticPressure(*this)
            );
        }

        //- Destructor
        virtual ~hydrostaticPressure();

    //- Update values
    virtual void update()
    {}

    //- Return pressure given shape function and nodes
    virtual scalar pressure
    (
        const labelList& nodeLabels,
        const scalarList& shape,
        const vector& n
    ) const;

    //- Return traction given shape function, nodes, and normal vector
    virtual vector traction
    (
        const labelList& nodeLabels,
        const scalarList& shape,
        const vector& n
    ) const;
};


} // End namespace boundaryTractions
} // End namespace Foam

#endif
