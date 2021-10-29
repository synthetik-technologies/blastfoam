/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryObjectLinearDamper.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBoundaryObjectRestraints
{
    defineTypeNameAndDebug(linearDamper, 0);

    addToRunTimeSelectionTable
    (
        immersedBoundaryObjectRestraint,
        linearDamper,
        motion
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectRestraints::linearDamper::linearDamper
(
    const polyMesh& mesh,
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    immersedBoundaryObjectRestraint(mesh, name, sDoFRBMRDict),
    coeff_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryObjectRestraints::linearDamper::~linearDamper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::immersedBoundaryObjectRestraints::linearDamper::restrain
(
    const immersedBoundaryObject& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintForce = -coeff_*motion.v();
    restraintMoment = Zero;

    if (motion.report())
    {
        Info<< " force " << restraintForce
            << endl;
    }
}


bool Foam::immersedBoundaryObjectRestraints::linearDamper::read
(
    const dictionary& sDoFRBMRDict
)
{
    immersedBoundaryObjectRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.lookup("coeff") >> coeff_;

    return true;
}


void Foam::immersedBoundaryObjectRestraints::linearDamper::write
(
    Ostream& os
) const
{
    writeEntry(os, "coeff", coeff_);
}


// ************************************************************************* //
