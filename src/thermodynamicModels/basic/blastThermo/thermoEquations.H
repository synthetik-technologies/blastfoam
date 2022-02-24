/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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
    Foam::fluidEEquation
    Foam::THEEquation

Description
    Classes used to solve for internal energy and temperature

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef thermoEquations_H
#define thermoEquations_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidBlastThermo.H"
#include "equation.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidEEquation Declaration
\*---------------------------------------------------------------------------*/

class fluidEEquation
:
    public equation
{
    //- Constant reference to the fluidBlastThermo
    const fluidBlastThermo& thermo_;

    //- Constant reference to the pressure
    const volScalarField& p_;

    //- Non-constant reference to the internal energy
    //  A pointer is used so it can be mutable
    mutable volScalarField* e_;

    //- Saved old internal energy
    scalar eOld_;

public:
    fluidEEquation(fluidBlastThermo& thermo)
    :
        equation(-great, great),
        thermo_(thermo),
        p_(thermo_.p()),
        e_(&thermo.he()),
        eOld_(0.0)
    {}

    //- Return the number of derivatives
    virtual label nDerivatives() const
    {
        return 1;
    }

    // Save the old energy
    void save(const label li)
    {
        eOld_ = (*e_)[li];
    }

    // Revert to the old energy
    void reset(const label li)
    {
        (*e_)[li] = eOld_;
    }

    //- Return the difference between the calulcated pressure and the current
    virtual scalar fx(const scalar e, const label li) const
    {
        (*e_)[li] = e;
        return thermo_.cellpRhoT(li, false) - p_[li];
    }

    //- Return the derivative of pressure w.r.t. energy
    virtual scalar dfdx(const scalar e, const label li) const
    {
        (*e_)[li] = e;
        return thermo_.celldpde(li);
    }
};


/*---------------------------------------------------------------------------*\
                        Class THEEquation Declaration
\*---------------------------------------------------------------------------*/


class THEEquation
:
    public equation
{
    //- Constant reference to the blastThermo
    const blastThermo& thermo_;

    //- Constant reference to the internal energy field
    const volScalarField& he_;

    //- Current patch index
    label patchi_;

public:
    THEEquation
    (
        fluidBlastThermo& thermo,
        const scalar TLow
    )
    :
        equation(-great, great),
        thermo_(thermo),
        he_(thermo.he()),
        patchi_(-1)
    {}

    //- Return the number of derivatives
    virtual label nDerivatives() const
    {
        return 1;
    }

    //- Set the current patch index
    void setPatch(const label patchi)
    {
        patchi_ = patchi;
    }

    //- Set the current patch index to -1
    void setCells()
    {
        patchi_ = -1;
    }

    //- Return the difference between the calculated internal energy and the
    //  current
    virtual scalar fx(const scalar T, const label li) const
    {
        return
            patchi_ == -1
          ? thermo_.cellHE(T, li) - he_[li]
          : thermo_.patchFaceHE(T, patchi_, li)
          - he_.boundaryField()[patchi_][li];
    }

    //- Return the derivative of temperature w.r.t. internal energy
    virtual scalar dfdx(const scalar T, const label li) const
    {
        return
            patchi_ == -1
          ? thermo_.cellCpv(T, li)
          : thermo_.patchFaceCpv(T, patchi_, li);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
