/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Find the reference cells per region. Regions given as cellSet names.

SourceFiles
    findRefCell.C

SeeAlso
    Foam::regionSplit

\*---------------------------------------------------------------------------*/

#ifndef findRefCells_H
#define findRefCells_H

#include "label.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Set reference cells for multi-region domains. Returns per region:
//  - whether a reference is needed
//  - if so a reference cell (or -1)
//  - and a reference value
void setRefCells
(
    const volScalarField& field,
    const volScalarField& fieldRef,
    const dictionary& dict,

    const label refCelli,
    const scalar refValue,

    boolList& needReference,
    labelList& refCells,
    scalarField& refValues,
    const bool forceReference = false
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
