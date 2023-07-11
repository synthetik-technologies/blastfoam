/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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

#include "lookupTables2D.H"
#include "linearLeastSquares.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * Scalar Functions * * * * * * * * * * * * * * //
template<>
Foam::scalar Foam::lookupTable2D<Foam::scalar>::reverseLookupX
(
    const scalar& fin,
    const scalar y
) const
{
    return solver(0, fin, y).solveUni();
}


template<>
Foam::scalar Foam::lookupTable2D<Foam::scalar>::reverseLookupY
(
    const scalar& fin,
    const scalar x
) const
{
    return solver(1, fin, x).solveUni();
}


namespace Foam
{
typedef HashSet<labelPair, labelPair::Hash<>> labelPairHashSet;
void addSamples
(
    const label i, const label j,
    const label nx, const label ny,
    List2D<bool>& added,
    labelPairHashSet& next
)
{
    if (i > 0 && !added(i-1, j))
    {
        next.insert({i-1, j});
    }
    if (j > 0 && !added(i, j-1))
    {
        next.insert({i, j-1});
    }
    if (i < nx-1 && !added(i+1, j))
    {
        next.insert({i+1, j});
    }
    if (j < ny-1 && !added(i, j+1))
    {
        next.insert({i, j+1});
    }

    if ((i > 0 && j > 0) && !added(i-1, j-1))
    {
        next.insert({i-1, j-1});
    }
    if ((i > 0 && j < ny-1) && !added(i-1, j+1))
    {
        next.insert({i-1, j+1});
    }
    if ((i < nx-1 && j > 0) && !added(i+1, j-1))
    {
        next.insert({i+1, j-1});
    }
    if ((i > nx-1 && j < ny-1) && !added(i+1, j+1))
    {
        next.insert({i+1, j+1});
    }
}


List2D<scalar> leastSquaresFit
(
    const scalarList& xSamples,
    const scalarList& ySamples,
    const scalarField& FSamples,
    const label minSamples,
    lookupTable2D<scalar>& table
)
{
    const scalarList& xMod = table.xMod();
    const scalarList& yMod = table.yMod();

    const scalarList& x = table.x();
    const scalarList& y = table.y();

    const label nx = x.size();
    const label ny = y.size();
    Field<vector2D> XSamples(xSamples.size());
    Field<vector2D> XModSamples(xSamples.size());
    forAll(XModSamples, i)
    {
        XSamples[i][0] = xSamples[i];
        XSamples[i][1] = ySamples[i];
        XModSamples[i][0] = table.modX()(xSamples[i]);
        XModSamples[i][1] = table.modY()(ySamples[i]);
    }

    List2D<labelHashSet> rMap(nx, ny);
    forAll(XSamples, samplei)
    {
        table.updateIndex(XSamples[samplei][0], XSamples[samplei][1]);
        const label i = table.i();
        const label j = table.j();
        rMap(i, j).insert(samplei);
        rMap(i+1, j).insert(samplei);
        rMap(i+1, j+1).insert(samplei);
        rMap(i, j+1).insert(samplei);
    }

    List2D<labelHashSet> fullMap(nx, ny);
    List2D<bool> added(nx, ny, false);
    for (label i = 0; i < nx; i++)
    {
        for (label j = 0; j < ny; j++)
        {
            added.List<bool>::operator=(false);
            labelPairHashSet prev({{i, j}});
            labelPairHashSet next;
            labelHashSet& mapi = fullMap(i, j);

            while (mapi.size() < minSamples)
            {
                next.clear();
                forAllConstIter
                (
                    labelPairHashSet,
                    prev,
                    iter
                )
                {
                    label ii = iter.key().first();
                    label jj = iter.key().second();
                    mapi += rMap(ii, jj);
                    added(ii, jj) = true;
                    addSamples
                    (
                        ii, jj,
                        rMap.m(), rMap.n(),
                        added,
                        next
                    );
                }
                prev = next;
            }
        }
    }

    List2D<scalar> f(nx, ny);
    linearLeastSquares solver;
    for (label i = 0; i < nx; i++)
    {
        for (label j = 0; j < ny; j++)
        {
            labelList indices(fullMap(i, j).toc());
            List<vector2D> subX(XModSamples, indices);
            Field<scalar> subF(FSamples, indices);

            autoPtr<scalarUnivariateEquation> eqn(solver.createEquation(subX, subF));
            f(i, j) = eqn->fX({xMod[i], yMod[j]}, 0);
        }
    }
    return f;
}
}

template<>
void Foam::lookupTable2D<Foam::scalar>::read
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name,
    const bool canRead
)
{
    xName_ = xName;
    yName_ = yName;
    const word scheme
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linearClamp")
    );

    scalarField x;
    {
        const dictionary& xDict = readComponent
        (
            dict,
            xName,
            modX_,
            x,
            canRead
        );
        setX(x, true);
        xInterpolator_ = interpolationWeight1D::New
        (
            xDict.found("interpolationScheme")
          ? xDict.lookup<word>("interpolationScheme")
          : dict.lookupOrDefault<word>
            (
                xName + "InterpolationScheme",
                scheme
            ),
            xModValues_,
            canRead
        );
        if (canRead)
        {
            xInterpolator_->validate();
        }
    }

    {
        scalarField y;
        const dictionary& yDict = readComponent
        (
            dict,
            yName,
            modY_,
            y,
            canRead
        );
        setY(y, true);
        yInterpolator_ = interpolationWeight1D::New
        (
            yDict.found("interpolationScheme")
          ? yDict.lookup<word>("interpolationScheme")
          : dict.lookupOrDefault<word>
            (
                yName + "InterpolationScheme",
                scheme
            ),
            yModValues_,
            canRead
        );
        if (canRead)
        {
            yInterpolator_->validate();
        }
    }

    List2D<scalar> data(xModValues_.size(), yModValues_.size());

    bool useLeastSquares = dict.lookupOrDefault("leastSquaresFit", false);
    if (dict.found(name) || dict.found(name + "File"))
    {
        if (!useLeastSquares && canRead)
        {
            if (dict.found(name))
            {
                dict.readIfPresent(name, data);
            }
            else
            {
                read2DTable
                (
                    dict.lookup<fileName>(name + "File"),
                    readDelim(dict, name + "Delim"),
                    data,
                    dict.lookupOrDefault<bool>(name + "FlipTable", false),
                    !canRead
                );
            }
        }
        mod_ = Modifier<scalar>::New
        (
            dict.lookupOrDefault<word>(name + "Mod", "none"),
            dict
        );
        mod_->readReal(dict, name + "IsReal");
    }
    else if (dict.isDict(name + "Coeffs"))
    {
        const dictionary& fDict(dict.subDict(name + "Coeffs"));
        mod_ = Modifier<scalar>::New
        (
            fDict.lookupOrDefault<word>("mod", "none"),
            fDict
        );
        mod_->readReal(fDict, "isReal");

        if (useLeastSquares || !canRead)
        {}
        else if (fDict.found(name))
        {
            fDict.readIfPresent(name, data);
        }
        else if (fDict.found("file"))
        {
            read2DTable
            (
                fDict.lookup<fileName>("file"),
                readDelim(fDict),
                data,
                fDict.lookupOrDefault<bool>("flipTable", false),
                !canRead
            );
        }
        else
        {
            FatalIOErrorInFunction(fDict)
                << "Neither the entry \"" << name << "\", "
                << " or a file was provided for construction" << endl
                << abort(FatalIOError);
        }
    }
    else if (!useLeastSquares || !canRead)
    {
        FatalIOErrorInFunction(dict)
            << "Neither the entry \"" << name << "\", "
            << " or the \""
            << name << "Coeffs\" subDictionary was found" << endl
            << abort(FatalIOError);
    }

    if (!mod_.valid())
    {
        mod_ = Modifier<scalar>::New
        (
            dict.lookupOrDefault<word>(name + "Mod", "none"),
            dict
        );
        mod_->readReal(dict, name + "IsReal");
    }

    if (useLeastSquares)
    {
        const dictionary& lsDict(dict.subDict("leastSquaresCoeffs"));
        scalarField sparseF;
        scalarField sparseX;
        scalarField sparseY;
        autoPtr<Modifier<scalar>> sparseXMod;
        readComponent
        (
            lsDict,
            xName,
            sparseXMod,
            sparseX,
            true
        );

        autoPtr<Modifier<scalar>> sparseYMod;
        readComponent
        (
            lsDict,
            yName,
            sparseYMod,
            sparseY,
            true
        );

        autoPtr<Modifier<scalar>> sparseFMod;
        readComponent
        (
            lsDict,
            name,
            sparseFMod,
            sparseF,
            true
        );
        mod_->Mod(sparseF);

        label minSamples = lsDict.lookupOrDefault<label>("minSamples", 4);
        data = leastSquaresFit
        (
            sparseX,
            sparseY,
            sparseF,
            minSamples,
            *this
        );
        mod_->Inv(data);
        mod_->setReal();
    }

    if
    (
        data.m() != xModValues_.size()
     || data.n() != yModValues_.size()
    )
    {
        FatalIOErrorInFunction(dict)
            << "Incompatible dimensions for table" << nl
            << "table size: "
            << data.m() << " x " << data.n() << nl
            << "x and y size: "
            << xModValues_.size() << " x " << yModValues_.size() << nl
            << abort(FatalIOError);
    }
    if (!mod_->isReal())
    {
        mod_->Inv(data);
        mod_->setReal();
    }
    setData(data, true);

    if (dict.found("rootSolver"))
    {
        this->solver
        (
            dict.lookup<word>("rootSolver"),
            dict
        );
    }
}

// ************************************************************************* //
