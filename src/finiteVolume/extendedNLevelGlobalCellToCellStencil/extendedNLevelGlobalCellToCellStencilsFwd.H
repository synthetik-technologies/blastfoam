/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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
    Foam::extendedNLevelGlobalCellToCellStencil

Description
    Creates an n-level cell stencil using global indicies

SourceFiles
    extendedNLevelGlobalCellToCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedNLevelGlobalCellToCellStencilsFwd_H
#define extendedNLevelGlobalCellToCellStencilsFwd_H

#include "CFCCellToCellStencil.H"
#include "CECCellToCellStencil.H"
#include "CPCCellToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class StencilType>
class extendedNLevelGlobalCellToCellStencil;

typedef extendedNLevelGlobalCellToCellStencil<CFCCellToCellStencil>
    extendedNLevelCFCCellToCellStencil;
typedef extendedNLevelGlobalCellToCellStencil<CECCellToCellStencil>
    extendedNLevelCECCellToCellStencil;
typedef extendedNLevelGlobalCellToCellStencil<CPCCellToCellStencil>
    extendedNLevelCPCCellToCellStencil;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
