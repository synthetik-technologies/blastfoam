/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Virtual base class for cohesive law.

\*---------------------------------------------------------------------------*/

#include "simpleCohesiveZoneLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleCohesiveZoneLaw, 0);
    defineRunTimeSelectionTable(simpleCohesiveZoneLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::simpleCohesiveZoneLaw> Foam::simpleCohesiveZoneLaw::New
(
    const word& simpleCohesiveZoneLawName,
    const dictionary& dict
)
{
    if (debug)
    {
        Info << "Selecting cohesive law: " << simpleCohesiveZoneLawName << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(simpleCohesiveZoneLawName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "simpleCohesiveZoneLaw::New(const word& simpleCohesiveZoneLawName, "
            "const dictionary& dict)"
        )   << "Unknown cohesive law " << simpleCohesiveZoneLawName
            << endl << endl
            << "Valid cohesive laws are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<simpleCohesiveZoneLaw>
    (
        cstrIter()(simpleCohesiveZoneLawName, dict)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::simpleCohesiveZoneLaw::simpleCohesiveZoneLaw
(
    const word& simpleCohesiveZoneLawName,
    const dictionary& dict
)
:
    simpleCohesiveZoneLawCoeffs_
    (
        dict.subDict(simpleCohesiveZoneLawName + "Coeffs")
    ),
    GIc_(simpleCohesiveZoneLawCoeffs_.lookup("GIc")),
    sigmaMax_(simpleCohesiveZoneLawCoeffs_.lookup("sigmaMax"))
{}


Foam::simpleCohesiveZoneLaw::simpleCohesiveZoneLaw
(
    const simpleCohesiveZoneLaw& cl
)
:
    refCount(),
    simpleCohesiveZoneLawCoeffs_(cl.simpleCohesiveZoneLawCoeffs_),
    GIc_(cl.GIc_),
    sigmaMax_(cl.sigmaMax_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleCohesiveZoneLaw::~simpleCohesiveZoneLaw()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleCohesiveZoneLaw::writeDict(Ostream& os) const
{
    os.writeKeyword(word(type() + "Coeffs"))
        << simpleCohesiveZoneLawCoeffs();
}


// ************************************************************************* //
