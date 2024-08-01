/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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
    Foam::boundaryTractions::idealGasPressure

Description
    Base class for constraining nodal displacement and velocity

SourceFiles
    displacementConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef idealGasPressureBoundaryTraction_H
#define idealGasPressureBoundaryTraction_H

#include "boundaryTraction.H"

namespace Foam
{
namespace boundaryTractions
{

class idealGasPressure
:
    public boundaryTraction
{
protected:

    //- Initial pressure
    scalar p0_;

    //- Reference volume
    scalar VRef_;

    //- Reference pressure
    scalar pRef_;

    //- Current pressure
    scalar p_;

    //- Calculate volume internal to a set of faces
    static scalar calcV
    (
        const pointField& pts,
        const faceList& faces
    );

    //- Copy constructor
    idealGasPressure(const idealGasPressure& igpbt);

public:

    //- Runtime type information
    TypeName("idealGasPressure")

        // Constructor
        idealGasPressure
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
                new idealGasPressure(*this)
            );
        }

        //- Destructor
        virtual ~idealGasPressure();

    //- Update values
    virtual void update();

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
