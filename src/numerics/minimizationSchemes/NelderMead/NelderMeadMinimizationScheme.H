/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::NelderMeadMinimizationScheme

Description
    Gradient free minimization scheme of Nelder-Mead

References:
    \verbatim
        "Algorithms for Optimization",
        Kochenderfer, M. J., Wheeler, T. A.,
        The MIT Press, Cambridge, Massachusetts, 2019.
    \endverbatim

SourceFiles
    NelderMeadMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef NelderMeadMinimizationScheme_H
#define NelderMeadMinimizationScheme_H

#include "minimizationScheme.H"
#include "typeInfo.H"
#include "Random.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class NelderMeadMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class NelderMeadMinimizationScheme
:
    public minimizationScheme
{
private:

    //- Coefficients
    scalar reflectionCoeff_;
    scalar expansionCoeff_;
    scalar contractionCoeff_;

    //- Return the mean of a list of points
    template<template <class> class ListListType, class ListType>
    void mean
    (
        const ListListType<ListType>& lst,
        ListType& avg,
        const label n
    ) const
    {
        if (n < 1)
        {
            return;
        }
        avg = lst[0];
        for (label i = 1; i < n; i++)
        {
            avg += lst[i];
        }
        avg /= scalar(n);
    }

    //- Return the mean of a list of points
    template<template <class> class ListListType, class ListType>
    void mean(const ListListType<ListType>& lst, ListType& avg) const
    {
        mean(lst, avg, lst.size());
    }


    template<template <class> class ListListType, class ListType>
    void variance
    (
        const ListListType<ListType>& lst,
        const ListType& avg,
        ListType& var
    ) const
    {
        var = sqr(lst[0] - avg);
        for (label i = 1; i < lst.size(); i++)
        {
            var += sqr(lst[i] - avg);
        }
        var /= scalar(lst.size());
    }


public:

    //- Runtime type information
    TypeName("NelderMead");

    // Constructors

        //- Construct for given minimizationScheme
        NelderMeadMinimizationScheme
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

    //- Destructor
    virtual ~NelderMeadMinimizationScheme()
    {}


    // Member Functions

        //- Solve the root
        virtual tmp<scalarField> minimize
        (
            const scalarList& x0,
            const scalarList& x1,
            const scalarList& x2,
            const label li
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
