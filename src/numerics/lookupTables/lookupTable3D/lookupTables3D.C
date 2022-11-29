/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "lookupTables3D.H"
#include "linearLeastSquares.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
typedef HashSet<FixedList<label, 3>, FixedList<label, 3>::Hash<>> labelTripHashSet;
void addSamples
(
    const label i, const label j, const label k,
    const label nx, const label ny, const label nz,
    List3D<bool>& added,
    labelTripHashSet& next
)
{
    if (i > 0 && !added(i-1, j, k))
    {
        next.insert({i-1, j, k});
    }
    if (j > 0 && !added(i, j-1, k))
    {
        next.insert({i, j-1, k});
    }
    if (k > 0 && !added(i, j, k-1))
    {
        next.insert({i, j, k-1});
    }
    if (i < nx-1 && !added(i+1, j, k))
    {
        next.insert({i+1, j, k});
    }
    if (j < ny-1 && !added(i, j+1, k))
    {
        next.insert({i, j+1, k});
    }
    if (k < nz-1 && !added(i, j, k+1))
    {
        next.insert({i, j, k+1});
    }

    if ((i > 0 && j > 0) && !added(i-1, j-1, k))
    {
        next.insert({i-1, j-1, k});
    }
    if ((i > 0 && k > 0) && !added(i-1, j, k-1))
    {
        next.insert({i-1, j, k-1});
    }
    if ((i > 0 && j < ny-1) && !added(i-1, j+1, k))
    {
        next.insert({i-1, j+1, k});
    }
    if ((i > 0 && k < nz-1) && !added(i-1, j, k+1))
    {
        next.insert({i-1, j, k+1});
    }
    if ((j > 0 && k > 0) && !added(i, j-1, k-1))
    {
        next.insert({i, j-1, k-1});
    }
    if ((j > 0 && k < nz-1) && !added(i, j-1, k+1))
    {
        next.insert({i, j-1, k+1});
    }

    if ((i < nx-1 && j > 0) && !added(i+1, j-1, k))
    {
        next.insert({i+1, j-1, k});
    }
    if ((i < nx-1 && k > 0) && !added(i+1, j, k-1))
    {
        next.insert({i+1, j, k-1});
    }
    if ((i > nx-1 && j < ny-1) && !added(i+1, j+1, k))
    {
        next.insert({i+1, j+1, k});
    }
    if ((i > nx-1 && k < nz-1) && !added(i+1, j, k+1))
    {
        next.insert({i+1, j, k+1});
    }
    if ((j < ny-1 && k > 0) && !added(i, j+1, k-1))
    {
        next.insert({i, j+1, k-1});
    }
    if ((j > ny-1 && k < nz-1) && !added(i, j+1, k+1))
    {
        next.insert({i, j+1, k+1});
    }

    if ((i > 0 && j > 0 && k > 0) && !added(i-1, j-1, k-1))
    {
        next.insert({i-1, j-1, k-1});
    }
    if ((i < nx-1 && j > 0 && k > 0) && !added(i+1, j-1, k-1))
    {
        next.insert({i+1, j-1, k-1});
    }
    if ((i > 0 && j < ny-1 && k > 0) && !added(i-1, j+1, k-1))
    {
        next.insert({i-1, j+1, k-1});
    }
    if ((i < nx-1 && j < ny-1 && k > 0) && !added(i+1, j+1, k-1))
    {
        next.insert({i+1, j+1, k-1});
    }

    if ((i > 0 && j > 0 && k < nz-1) && !added(i-1, j-1, k+1))
    {
        next.insert({i-1, j-1, k+1});
    }
    if ((i < nx-1 && j > 0 && k < nz-1) && !added(i+1, j-1, k+1))
    {
        next.insert({i+1, j-1, k+1});
    }
    if ((i > 0 && j < ny-1 && k < nz-1) && !added(i-1, j+1, k+1))
    {
        next.insert({i-1, j+1, k+1});
    }
    if ((i < nx-1 && j > ny-1 && k < nz-1) && !added(i+1, j+1, k+1))
    {
        next.insert({i+1, j+1, k+1});
    }
}

List3D<scalar> leastSquaresFit
(
    const scalarList& xSamples,
    const scalarList& ySamples,
    const scalarList& zSamples,
    const scalarField& FSamples,
    const label minSamples,
    lookupTable3D<scalar>& table
)
{
    const scalarList& xMod = table.xMod();
    const scalarList& yMod = table.yMod();
    const scalarList& zMod = table.zMod();

    const scalarList& x = table.x();
    const scalarList& y = table.y();
    const scalarList& z = table.z();

    const label nx = x.size();
    const label ny = y.size();
    const label nz = z.size();

    Field<point> XSamples(xSamples.size());
    Field<point> XModSamples(xSamples.size());
    forAll(XModSamples, i)
    {
        XSamples[i][0] = xSamples[i];
        XSamples[i][1] = ySamples[i];
        XSamples[i][2] = zSamples[i];
        XModSamples[i][0] = table.modX()(xSamples[i]);
        XModSamples[i][1] = table.modY()(ySamples[i]);
        XModSamples[i][2] = table.modZ()(zSamples[i]);
    }

    // List3D<labelHashSet> rMap(nx, ny, nz);
    // forAll(XSamples, samplei)
    // {
    //     table.updateIndex
    //     (
    //         XSamples[samplei][0],
    //         XSamples[samplei][1],
    //         XSamples[samplei][2]
    //     );
    //     const label i = table.i();
    //     const label j = table.j();
    //     const label k = table.k();
    //     rMap(i, j, k).insert(samplei);
    //     rMap(i+1, j, k).insert(samplei);
    //     rMap(i+1, j+1, k).insert(samplei);
    //     rMap(i, j+1, k).insert(samplei);
    //     rMap(i, j, k+1).insert(samplei);
    //     rMap(i+1, j, k+1).insert(samplei);
    //     rMap(i+1, j+1, k+1).insert(samplei);
    //     rMap(i, j+1, k+1).insert(samplei);
    // }
    //
    // List3D<labelHashSet> fullMap(nx, ny, nz);
    // List3D<bool> added(nx, ny, nz, false);
    // for (label i = 0; i < nx; i++)
    // {
    //     for (label j = 0; j < ny; j++)
    //     {
    //         for (label k = 0; k < nz; k++)
    //         {
    //             added.List<bool>::operator=(false);
    //             labelTripHashSet prev({{i, j, k}});
    //             labelTripHashSet next;
    //             labelHashSet& mapi = fullMap(i, j, k);
    //
    //             while (mapi.size() < minSamples)
    //             {
    //                 next.clear();
    //                 forAllConstIter
    //                 (
    //                     labelTripHashSet,
    //                     prev,
    //                     iter
    //                 )
    //                 {
    //                     label ii = iter.key()[0];
    //                     label jj = iter.key()[1];
    //                     label kk = iter.key()[2];
    //                     mapi += rMap(ii, jj, kk);
    //                     added(ii, jj, kk) = true;
    //                     addSamples
    //                     (
    //                         ii, jj, kk,
    //                         rMap.m(), rMap.n(), rMap.l(),
    //                         added,
    //                         next
    //                     );
    //                 }
    //                 prev = next;
    //             }
    //         }
    //     }
    // }

    List3D<labelHashSet> fullMap(nx, ny, nz);
    treeBoundBox overallBb(XModSamples);
    indexedOctree<treeDataPoint> tree
    (
        treeDataPoint(XModSamples),
        overallBb,
        10,
        100,
        100.0
    );

    treeBoundBox bb(overallBb);
    for (label i = 0; i < nx; i++)
    {
        overallBb.min().x() = i > 0 ? xMod[i-1] : xMod[i];
        overallBb.max().x() = i < nx-1 ? xMod[i+1] : xMod[i];

        for (label j = 0; j < ny; j++)
        {
            overallBb.min().y() = j > 0 ? yMod[j-1] : yMod[j];
            overallBb.max().y() = j < ny-1 ? yMod[j+1] : yMod[j];
            for (label k = 0; k < nz; k++)
            {
                overallBb.min().z() = k > 0 ? zMod[k-1] : zMod[k];
                overallBb.max().z() = k < nz-1 ? zMod[k+1] : zMod[k];
                bb = overallBb;

                labelList indices(tree.findBox(bb));
                while (indices.size() < minSamples)
                {
                    bb.min() -= bb.span()/2.0;
                    bb.max() += bb.span()/2.0;
                    indices = move(tree.findBox(bb));
                }
                fullMap(i, j, k).insert(indices);
            }
        }
    }


    List3D<scalar> f(nx, ny, nz);
    linearLeastSquares solver;
    for (label i = 0; i < nx; i++)
    {
        for (label j = 0; j < ny; j++)
        {
            for (label k = 0; k < nz; k++)
            {
                labelList indices(fullMap(i, j, k).toc());
                List<vector> subX(XModSamples, indices);
                Field<scalar> subF(FSamples, indices);
                Field<scalar> weights(subF.size());
                forAll(weights, samplei)
                {
                    weights[samplei] =
                        1.0
                       /max
                        (
                            mag
                            (
                                subX[samplei]
                              - vector(xMod[i], yMod[j], zMod[k])
                            ),
                            small
                        );
                }

                autoPtr<scalarUnivariateEquation> eqn
                (
                    solver.createEquation(subX, subF, weights)
                );
                f(i, j, k) = eqn->fX({xMod[i], yMod[j], zMod[k]}, 0);
            }
        }
    }
    return f;
}
}


template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupX
(
    const scalar& fin,
    const scalar y,
    const scalar z
) const
{
    return solver(0, fin, y, z).solveUni();
}


template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupY
(
    const scalar& fin,
    const scalar x,
    const scalar z
) const
{
    return solver(1, fin, x, z).solveUni();
}


template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupZ
(
    const scalar& fin,
    const scalar x,
    const scalar y
) const
{

    return solver(2, fin, x, y).solveUni();
}


template<>
void Foam::lookupTable3D<Foam::scalar>::read
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& zName,
    const word& name,
    const bool canRead
)
{
    xName_ = xName;
    yName_ = yName;
    zName_ = zName;
    fName_ = name;

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
        xInterpolator_->validate();
    }

    scalarField y;
    {
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
        yInterpolator_->validate();
    }

    scalarField z;
    {
        const dictionary& zDict = readComponent
        (
            dict,
            zName,
            modZ_,
            z,
            canRead
        );
        setZ(z, true);
        zInterpolator_ = interpolationWeight1D::New
        (
            zDict.found("interpolationScheme")
          ? zDict.lookup<word>("interpolationScheme")
          : dict.lookupOrDefault<word>
            (
                zName + "InterpolationScheme",
                scheme
            ),
            zModValues_,
            canRead
        );
        zInterpolator_->validate();
    }

    List3D<scalar> data
    (
        xModValues_.size(),
        yModValues_.size(),
        zModValues_.size()
    );

    bool useLeastSquares = dict.lookupOrDefault("leastSquaresFit", false);
    if (dict.found(name))
    {
        if (!useLeastSquares)
        {
            dict.readIfPresent(name, data);
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

        if (useLeastSquares)
        {}
        else if (fDict.found(name))
        {
            fDict.readIfPresent(name, data);
        }
        else if (fDict.found("file"))
        {
            read3DTable
            (
                fDict.lookup<fileName>("file"),
                dict.lookupOrDefault<char>("delim", ','),
                dict.lookupOrDefault<char>("rowDelim", ';'),
                data,
                dict.lookupOrDefault<Switch>("flipTable", true),
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
    else if (!useLeastSquares)
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
        scalarField sparseZ;
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

        autoPtr<Modifier<scalar>> sparseZMod;
        readComponent
        (
            lsDict,
            zName,
            sparseZMod,
            sparseZ,
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
            sparseZ,
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
     || data.l() != zModValues_.size()
    )
    {
        FatalIOErrorInFunction(dict)
            << "Incompatible dimensions for table" << nl
            << "table size: "
            << data.m() << " x "
            << data.n() << " x "
            << data.l() << nl
            << "x and y size: "
            << xModValues_.size() << " x "
            << yModValues_.size() << " z "
            << yModValues_.size() << nl
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
