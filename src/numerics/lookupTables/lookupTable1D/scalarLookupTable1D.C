/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
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

#include "scalarLookupTable1D.H"
#include "tableReader.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarLookupTable1D::scalarLookupTable1D()
{}


Foam::scalarLookupTable1D::scalarLookupTable1D
(
    const dictionary& dict,
    const word& xName,
    const word& name
)
{
    read(dict, xName, name);
}


Foam::scalarLookupTable1D::scalarLookupTable1D
(
    const Field<scalar>& x,
    const Field<scalar>& data,
    const word& mod,
    const word& xMod,
    const word& interpolationScheme,
    const bool isReal
)
{
    set(x, data, mod, xMod, interpolationScheme, isReal);
}


Foam::scalarLookupTable1D::scalarLookupTable1D
(
    const Field<scalar>& x,
    const word& xMod,
    const word& interpolationScheme,
    const bool isReal
)
{
    setX(x, xMod, interpolationScheme, isReal);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::scalarLookupTable1D::~scalarLookupTable1D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::scalarLookupTable1D::reverseLookup(const scalar yin) const
{
#ifdef FULL_DEBUG
    if (!modFunc_)
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    scalar y(modFunc_(yin));
    if (y < data_[0])
    {
        index_ = 0;
    }
    if (y > data_.last())
    {
        index_ = data_.size() - 2;
    }
    else
    {
        for (index_ = 0; index_ < data_.size(); index_++)
        {
            if (y < data_[index_])
            {
                index_--;
                break;
            }
        }
    }

    const scalar& ym(data_[index_]);
    const scalar& yp(data_[index_+1]);

    f_ = linearWeight(y, ym, yp);

    return invModXFunc_
    (
        xModValues_[index_]
      + f_*(xModValues_[index_+1] - xModValues_[index_])
    );
}


void Foam::scalarLookupTable1D::read
(
    const dictionary& dict,
    const word& xName,
    const word& name
)
{
    Info<<"read"<<endl;
    if (dict.found("file2D"))
    {
        word interpolationScheme
        (
            dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
        );
        Field<Field<scalar>> tdata;
        read2DTable
        (
            dict.lookup<fileName>("file2D"),
            dict.lookupOrDefault<string>("delim", ","),
            tdata,
            dict.lookupOrDefault<Switch>("flipTable", true),
            true
        );

        Field<scalar> xValues(tdata.size());
        Field<scalar> yValues(tdata.size());
        label xi = dict.lookupOrDefault<label>(xName + "Col", 0);
        label yi = dict.lookupOrDefault<label>(name + "Col", 1);
        forAll(tdata, coli)
        {
            xValues[coli] = tdata[coli][xi];
            yValues[coli] = tdata[coli][yi];
        }

        Switch isReal(dict.lookupOrDefault<Switch>("isReal", true));
        word mod(dict.lookupOrDefault<word>(name + "Mod", "none"));
        word xMod(dict.lookupOrDefault<word>(xName + "Mod", "none"));

        set
        (
            xValues,
            yValues,
            mod,
            xMod,
            interpolationScheme,
            isReal
        );
        return;
    }

    lookupTable1D<scalar>::read(dict, xName, name);
}


// ************************************************************************* //
