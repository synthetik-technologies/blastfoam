/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

inline const Foam::Time& Foam::timeIntegrator::time() const
{
    return obr_.time();
}


inline const Foam::timeIntegratorCoeffs& Foam::timeIntegrator::coeffs() const
{
    return coeffs_();
}



inline bool Foam::timeIntegrator::restart() const
{
    return restart_;
}


inline const Foam::List<Foam::scalarList>&
Foam::timeIntegrator::as() const
{
    return as_;
}


inline const Foam::List<Foam::scalarList>&
Foam::timeIntegrator::bs() const
{
    return bs_;
}


inline const Foam::scalarList& Foam::timeIntegrator::a() const
{
    return as_[stepi_];
}


inline const Foam::scalarList& Foam::timeIntegrator::b() const
{
    return bs_[stepi_];
}


inline Foam::label Foam::timeIntegrator::nSteps() const
{
    return as_.size();
}


inline Foam::label Foam::timeIntegrator::nOld() const
{
    return nOld_;
}


inline Foam::label Foam::timeIntegrator::getOldIndex(const label stepi) const
{
    return oldIs_[stepi];
}


inline Foam::label Foam::timeIntegrator::nDelta() const
{
    return nDelta_;
}


inline Foam::label Foam::timeIntegrator::getDeltaIndex(const label stepi) const
{
    return deltaIs_[stepi];
}


inline Foam::label Foam::timeIntegrator::step() const
{
    return stepi_;
}


inline Foam::scalar Foam::timeIntegrator::f() const
{
    return f_[stepi_];
}


inline Foam::scalar Foam::timeIntegrator::f0() const
{
    return f0_[stepi_];
}


inline bool Foam::timeIntegrator::firstStep() const
{
    return stepi_ == 0;
}


inline bool Foam::timeIntegrator::finalStep() const
{
    return stepi_ == coeffs_->nSteps()-1;
}


inline Foam::labelList Foam::timeIntegrator::oldSaveMap() const
{
    return coeffs_->oldSaveMap();
}


inline Foam::labelList Foam::timeIntegrator::deltaSaveMap() const
{
    return coeffs_->deltaSaveMap();
}


inline Foam::scalar Foam::timeIntegrator::maxCo() const
{
    return coeffs_->maxCo();
}


// ************************************************************************* //
