/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 Alberto Passalacqua
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "PDFODETransportModel.H"


// * * * * * * * * * * * * Public Member functions * * * * * * * * * * * * * //

template<class MomentType>
void Foam::PDFODETransportModel::saveFields
(
    const label stepi,
    const MomentType& newFields,
    PtrList<PtrList<volScalarField>>& storage
)
{
    forAll(newFields, mi)
    {
        storage[stepi].set
        (
            mi,
            new volScalarField(newFields[mi])
        );
    }
}

// ************************************************************************* //
