/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#ifdef OPENFOAMFOUNDATION

template<class SourcePatch, class TargetPatch>
inline Foam::label
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::singlePatchProc() const
{
    return singlePatchProc_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalar
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::
lowWeightCorrection() const
{
    return lowWeightCorrection_;
}


template<class SourcePatch, class TargetPatch>
inline bool
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::
applyLowWeightCorrection() const
{
    return lowWeightCorrection_ > 0;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::srcMagSf() const
{
    return srcMagSf_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::labelListList&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::srcAddress() const
{
    return srcAddress_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarListList&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::srcWeights() const
{
    return srcWeights_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarListList&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::srcWeights()
{
    return srcWeights_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::srcWeightsSum() const
{
    return srcWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarField&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::srcWeightsSum()
{
    return srcWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::mapDistribute&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::srcMap() const
{
    return srcMapPtr_();
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::tgtMagSf() const
{
    return tgtMagSf_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::labelListList&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::tgtAddress() const
{
    return tgtAddress_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarListList&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::tgtWeights() const
{
    return tgtWeights_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarListList&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::tgtWeights()
{
    return tgtWeights_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::tgtWeightsSum() const
{
    return tgtWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarField&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::tgtWeightsSum()
{
    return tgtWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::mapDistribute&
Foam::newAMIInterpolation<SourcePatch, TargetPatch>::tgtMap() const
{
    return tgtMapPtr_();
}


#endif // end of #ifdef OPENFOAMFOUNDATION

// ************************************************************************* //
