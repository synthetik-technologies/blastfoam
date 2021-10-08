/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License

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

inline const Foam::primitivePatchInterpolation&
Foam::mappedPatchSelector::pointInterpolator() const
{
    if (!pointInterpolatorPtr_.valid())
    {
        pointInterpolatorPtr_.set
        (
            new primitivePatchInterpolation(patch_)
        );
    }

    return pointInterpolatorPtr_();
}


inline const Foam::polyMesh& Foam::mappedPatchSelector::mesh() const
{
    return patch_.boundaryMesh().mesh();
}


inline const Foam::polyPatch& Foam::mappedPatchSelector::patch() const
{
    return patch_;
}


// ************************************************************************* //
