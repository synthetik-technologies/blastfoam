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

InClass
    mechanicalLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<mechanicalLaw> mechanicalLaw::NewLinGeomMechLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
{
    const word mechTypeName(dict.lookup("type"));

    Info<< "Selecting mechanical law " << mechTypeName << endl;

    linGeomMechLawConstructorTable::iterator cstrIter =
        linGeomMechLawConstructorTablePtr_->find(mechTypeName);

    if (cstrIter == linGeomMechLawConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "mechanicalLaw::New(\n"
            "    const word& name,\n"
            "    const fvMehs& mesh,\n"
            "    const dictionary& dict,\n"
            "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
            ")",
            dict
        )   << "Unknown mechanicalLaw type "
            << mechTypeName << endl << endl
            << "Valid linearGeometry mechanicalLaws are : " << endl
            << linGeomMechLawConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<mechanicalLaw>(cstrIter()(name, mesh, dict, nonLinGeom));
}


autoPtr<mechanicalLaw> mechanicalLaw::NewNonLinGeomMechLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
{
    const word mechTypeName(dict.lookup("type"));

    Info<< "Selecting mechanical law " << mechTypeName << endl;

    nonLinGeomMechLawConstructorTable::iterator cstrIter =
        nonLinGeomMechLawConstructorTablePtr_->find(mechTypeName);

    if (cstrIter == nonLinGeomMechLawConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "mechanicalLaw::New(\n"
            "    const word& name,\n"
            "    const fvMehs& mesh,\n"
            "    const dictionary& dict,\n"
            "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
            ")",
            dict
        )   << "Unknown mechanicalLaw type "
            << mechTypeName << endl << endl
            << "Valid nonLinearGeometry mechanicalLaws are : " << endl
            << nonLinGeomMechLawConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<mechanicalLaw>(cstrIter()(name, mesh, dict, nonLinGeom));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
