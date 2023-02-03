/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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
    Foam::displacementRelaxations::QNBase

Description
    QNBase relaxation

SourceFiles
    QNBaseDisplacementRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef QNBaseDisplacementRelaxation_H
#define QNBaseDisplacementRelaxation_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "displacementRelaxation.H"

namespace Foam
{
namespace displacementRelaxations
{

/*---------------------------------------------------------------------------*\
                           Class QNBase Declaration
\*---------------------------------------------------------------------------*/

class QNBase
:
    public displacementRelaxation
{
    scalar initRelaxFactor_;
    label nOldTimes_;

    List<vectorField> residuals_;
    List<vectorField> prevResiduals_;
    List<vectorField> oldResiduals_;

    List<DynamicList<vectorField>> pointsV_;
    List<DynamicList<vectorField>> pointsW_;
    List<DynamicList<scalar>> times_;

    //- Update saved fields
    void updateFields(const label iter, pointVectorField&);

public:

    TypeName("QNBase");

    // Constructor
    QNBase(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~QNBase();


    // Member functions

        //- Update residuals
        virtual void updateError(const pointVectorField&);

        //- Relax the field
        virtual void relax(const label iter, pointVectorField&);

        //- Clear residuals
        virtual void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace displacementRelaxations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
