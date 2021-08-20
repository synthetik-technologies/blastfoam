/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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

#include "phasePair.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH
(
    const volScalarField& d
) const
{
    return
        mag(dispersed().rho() - continuous().rho())
       *mag(g())
       *sqr(d)
       /sigma();
}

Foam::scalar Foam::phasePair::cellEoH
(
    const label celli,
    const scalar& d
) const
{
    return
        mag(dispersed().rho()[celli] - continuous().rho()[celli])
       *mag(g().value())
       *sqr(d)
       /cellsigma(celli);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePair::phasePair
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const bool ordered
)
:
    phasePairKey(phase1.name(), phase2.name(), ordered),
    phase1_(phase1),
    phase2_(phase2),
    g_(phase1.mesh().lookupObject<uniformDimensionedVectorField>("g"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phasePair::~phasePair()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::phaseModel& Foam::phasePair::dispersed() const
{
    FatalErrorInFunction
        << "Requested dispersed phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


const Foam::phaseModel& Foam::phasePair::continuous() const
{
    FatalErrorInFunction
        << "Requested continuous phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


Foam::word Foam::phasePair::name() const
{
    word name2(phase2().name());
    name2[0] = toupper(name2[0]);
    return phase1().name() + "And" + name2;
}


Foam::word Foam::phasePair::otherName() const
{
    word name1(first());
    name1[0] = toupper(name1[0]);
    return second() + "And" + name1;
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::rho() const
{
    return phase1()*phase1().rho() + phase2()*phase2().rho();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::magUr
(
    const label nodei,
    const label nodej
) const
{
    return mag(phase1().U(nodei) - phase2().U(nodej));
}


Foam::tmp<Foam::volVectorField> Foam::phasePair::Ur
(
    const label nodei,
    const label nodej
) const
{
    return dispersed().U(nodei) - continuous().U(nodej);
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Re
(
    const label nodei, const label nodej
) const
{
    return magUr(nodei, nodej)*dispersed().d(nodei)/continuous().nu();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Pr
(
    const label nodei,
    const label nodej
) const
{
    Info<<min(continuous().kappa())<<" "<<endl;
    return
         continuous().nu()
        *continuous().Cv()
        *continuous().rho()
        /continuous().kappa();
}

Foam::tmp<Foam::volScalarField> Foam::phasePair::We
(
    const label nodei, const label nodej
) const
{
    return
        sqr(magUr(nodei, nodej))*dispersed().d(nodei)
        *continuous().rho()/sigma_;
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Eo
(
    const label nodei,
    const label nodej
) const
{
    return EoH(dispersed().d(nodei));
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH1
(
    const label nodei,
    const label nodej
) const
{
    return
        EoH
        (
            dispersed().d(nodei)
           *cbrt(1.0 + 0.163*pow(Eo(nodei, nodej), 0.757))
        );
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH2
(
    const label nodei,
    const label nodej
) const
{
    return EoH(dispersed().d(nodei)/cbrt(E(nodei, nodej)));
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::sigma() const
{
    NotImplemented;
    return phase1();
//         phase1().fluid().lookupSubModel<surfaceTensionModel>
//         (
//             phasePair(phase1(), phase2())
//         ).sigma();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Mo() const
{
    return
        mag(g())
       *continuous().nu()
       *pow3(continuous().nu()*continuous().rho()/sigma());
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Ta
(
    const label nodei,
    const label nodej
) const
{
    return Re(nodei, nodej)*pow(Mo(), 0.23);
}

Foam::tmp<Foam::volScalarField> Foam::phasePair::E
(
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Requested aspect ratio of the dispersed phase in an unordered pair"
        << exit(FatalError);

    return phase1();
}


Foam::scalar Foam::phasePair::cellrho(const label celli) const
{
    return
        phase1()[celli]*phase1().rho()[celli]
      + phase2()[celli]*phase2().rho()[celli];
}


Foam::scalar Foam::phasePair::cellmagUr
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return mag(phase1().U(nodei)[celli] - phase2().U(nodej)[celli]);
}


Foam::vector Foam::phasePair::cellUr
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return dispersed().U(nodei)[celli] - continuous().U(nodej)[celli];
}


Foam::scalar Foam::phasePair::cellRe
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return
        cellmagUr(celli, nodei, nodej)
       *dispersed().celld(celli, nodei)
       /continuous().cellnu(celli);
}


Foam::scalar Foam::phasePair::cellPr
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return
         continuous().cellnu(celli)
        *continuous().cellCv(celli)
        *continuous().rho()[celli]
        /continuous().cellkappa(celli);
}

Foam::scalar Foam::phasePair::cellWe
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return
        sqr(cellmagUr(celli, nodei, nodej))*dispersed().celld(celli, nodei)
        *continuous().rho()[celli]/sigma_.value();
}


Foam::scalar Foam::phasePair::cellEo
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return cellEoH(celli, dispersed().celld(celli, nodei));
}


Foam::scalar Foam::phasePair::cellEoH1
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return
        cellEoH
        (
            celli,
            dispersed().celld(celli, nodei)
           *cbrt(1.0 + 0.163*pow(cellEo(celli, nodei, nodej), 0.757))
        );
}


Foam::scalar Foam::phasePair::cellEoH2
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return
        cellEoH
        (
            celli,
            dispersed().celld(celli, nodei)/cbrt(cellE(celli, nodei, nodej))
        );
}


Foam::scalar Foam::phasePair::cellsigma(const label celli) const
{
    NotImplemented;
    return phase1()[celli];
//         phase1().fluid().lookupSubModel<surfaceTensionModel>
//         (
//             phasePair(phase1(), phase2())
//         ).sigma();
}


Foam::scalar Foam::phasePair::cellMo(const label celli) const
{
    return
        mag(g().value())
       *continuous().cellnu(celli)
       *pow3
        (
            continuous().cellnu(celli)
           *continuous().rho()[celli]
           /cellsigma(celli)
        );
}


Foam::scalar Foam::phasePair::cellTa
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    return cellRe(celli, nodei, nodej)*pow(cellMo(celli), 0.23);
}

Foam::scalar Foam::phasePair::cellE
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    FatalErrorInFunction
        << "Requested aspect ratio of the dispersed phase in an unordered pair"
        << exit(FatalError);

    return phase1()[celli];
}
// ************************************************************************* //
