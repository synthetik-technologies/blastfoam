/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "blastChemicalEquilibrium.H"
#include "thermodynamicConstants.H"
#include "SVD.H"
#include <cctype>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::List<Foam::Tuple2<Foam::word, Foam::scalar>>
Foam::blastChemicalEquilibrium<ThermoType>::parseFormula(const word& name)
{
    const std::string f(name.substr(0, name.find('(')));

    List<Tuple2<word, scalar>> comp;
    std::string::size_type i = 0;
    while (i < f.size())
    {
        if (!std::isupper(f[i]))
        {
            FatalErrorInFunction
                << "Cannot parse the elements of specie " << name << nl
                << "    Give its composition in the elements dictionary"
                << exit(FatalError);
        }

        std::string el(1, f[i++]);
        while (i < f.size() && std::islower(f[i]))
        {
            el += f[i++];
        }

        std::string::size_type j = i;
        while (j < f.size() && (std::isdigit(f[j]) || f[j] == '.'))
        {
            j++;
        }
        const scalar count = j > i ? std::atof(f.substr(i, j - i).c_str()) : 1;
        i = j;

        bool found = false;
        forAll(comp, ci)
        {
            if (comp[ci].first() == el)
            {
                comp[ci].second() += count;
                found = true;
            }
        }
        if (!found)
        {
            comp.append(Tuple2<word, scalar>(word(el), count));
        }
    }

    return comp;
}


template<class ThermoType>
void Foam::blastChemicalEquilibrium<ThermoType>::gibbs
(
    const scalar T,
    scalarField& g
) const
{
    using constant::thermodynamic::RR;

    forAll(thermos_, k)
    {
        g[k] = thermos_[k].W()*thermos_[k].Gstd(T)/(RR*T);
    }
}


template<class ThermoType>
void Foam::blastChemicalEquilibrium<ThermoType>::activeSets
(
    const scalarField& b,
    boolList& activeEl,
    boolList& active
) const
{
    const scalar bRef = sum(b);
    activeEl.setSize(elements_.size());
    forAll(b, i)
    {
        activeEl[i] = b[i] > 1e-20*bRef;
    }

    active.setSize(nSpecie_);
    forAll(a_, k)
    {
        active[k] = true;
        forAll(a_[k], i)
        {
            if (a_[k][i] > 0 && !activeEl[i])
            {
                active[k] = false;
            }
        }
    }
}


template<class ThermoType>
void Foam::blastChemicalEquilibrium<ThermoType>::gasAmounts
(
    const scalar lnc0,
    const scalarField& g,
    const boolList& active,
    const scalarField& lambda,
    scalarField& n
) const
{
    forAll(n, k)
    {
        if (active[k] && !condensed_[k])
        {
            scalar x = lnc0 - g[k];
            forAll(lambda, i)
            {
                x += a_[k][i]*lambda[i];
            }
            n[k] = exp(min(x, 600.0));
        }
        else
        {
            n[k] = 0;
        }
    }
}


template<class ThermoType>
void Foam::blastChemicalEquilibrium<ThermoType>::initialLambda
(
    const scalar T,
    const scalar V,
    const scalarField& n0,
    const scalarField& b,
    const boolList& activeEl,
    const boolList& active,
    scalarField& lambda
) const
{
    using constant::thermodynamic::RR;
    using constant::thermodynamic::Pstd;

    // Weighted least-squares fit of the element potentials to the gas
    // composition: a_k.lambda = ln(n_k/c0) + g_k
    scalarField g(nSpecie_);
    gibbs(T, g);
    const scalar lnc0 = log(V*Pstd/(RR*T));
    const scalar nFloor = 1e-10*sum(b);

    labelList el;
    forAll(activeEl, i)
    {
        if (activeEl[i])
        {
            el.append(i);
        }
    }

    scalarSquareMatrix M(el.size(), Zero);
    List<scalar> rhs(el.size(), 0.0);
    forAll(n0, k)
    {
        if (!active[k] || condensed_[k])
        {
            continue;
        }
        const scalar nk = max(n0[k], 0) + nFloor;
        const scalar t = log(nk) - lnc0 + g[k];
        forAll(el, p)
        {
            const scalar ap = a_[k][el[p]];
            rhs[p] += nk*ap*t;
            forAll(el, q)
            {
                M(p, q) += nk*ap*a_[k][el[q]];
            }
        }
    }

    // Small regularisation for elements carried by few gas species
    forAll(el, p)
    {
        M(p, p) *= 1 + 1e-12;
        M(p, p) += vSmall;
    }
    LUsolve(M, rhs);

    lambda = 0;
    forAll(el, p)
    {
        lambda[el[p]] = rhs[p];
    }
}


template<class ThermoType>
bool Foam::blastChemicalEquilibrium<ThermoType>::solveKKT
(
    const scalarField& n,
    const scalarField& b,
    const boolList& activeEl,
    const labelList& C,
    const scalarField& r,
    const scalarField& rc,
    scalarField& dLambda,
    scalarField& nC
) const
{
    labelList el;
    forAll(activeEl, i)
    {
        if (activeEl[i])
        {
            el.append(i);
        }
    }
    const label nA = el.size();
    const label N = nA + C.size();

    scalarSquareMatrix M(N, Zero);
    List<scalar> x(N, 0.0);

    forAll(n, k)
    {
        if (n[k] > 0 && !condensed_[k])
        {
            forAll(el, p)
            {
                const scalar nap = n[k]*a_[k][el[p]];
                if (nap > 0)
                {
                    forAll(el, q)
                    {
                        M(p, q) += nap*a_[k][el[q]];
                    }
                }
            }
        }
    }
    forAll(C, c)
    {
        forAll(el, p)
        {
            M(p, nA + c) = a_[C[c]][el[p]];
            M(nA + c, p) = a_[C[c]][el[p]];
        }
        x[nA + c] = rc[c];
    }
    forAll(el, p)
    {
        x[p] = r[el[p]];
    }

    // Scaling: element rows by the element abundance, condensed amounts by
    // the abundance of their limiting element, so that trace elements are
    // resolved as well as the major ones
    scalarList sC(C.size(), great);
    forAll(C, c)
    {
        forAll(el, p)
        {
            const scalar ap = a_[C[c]][el[p]];
            if (ap > 0)
            {
                sC[c] = min(sC[c], b[el[p]]/ap);
            }
        }
    }
    forAll(el, p)
    {
        const scalar s = 1.0/b[el[p]];
        for (label q = 0; q < N; q++)
        {
            M(p, q) *= s;
        }
        x[p] *= s;
    }
    forAll(C, c)
    {
        for (label p = 0; p < N; p++)
        {
            M(p, nA + c) *= sC[c];
        }
    }

    // A non-finite or overflowing system (e.g. from an extreme cell state)
    // fails the solve, so that the cell is left frozen, instead of trapping
    // a floating-point exception in the SVD
    for (label p = 0; p < N; p++)
    {
        if (!std::isfinite(x[p]) || mag(x[p]) > 1e100)
        {
            return false;
        }
        for (label q = 0; q < N; q++)
        {
            if (!std::isfinite(M(p, q)) || mag(M(p, q)) > 1e100)
            {
                return false;
            }
        }
    }

    // Pseudo-inverse: directions not determined by significant species
    // (e.g. the potentials of elements held in a condensed phase, whose gas
    // species are negligible) are given no step. In the scaled system a
    // determined direction has a singular value of order one, so the cut-off
    // is absolute (a cut-off relative to the largest singular value would
    // also remove determined directions when the gas amounts of an element
    // overshoot its abundance during the iterations)
    const SVD svd(M);
    const scalarRectangularMatrix& U = svd.U();
    const scalarRectangularMatrix& V = svd.V();
    const scalarDiagonalMatrix& S = svd.S();
    List<scalar> y(N, 0.0);
    for (label j = 0; j < N; j++)
    {
        if (S[j] > 1e-12)
        {
            scalar uTx = 0;
            for (label q = 0; q < N; q++)
            {
                uTx += U(q, j)*x[q];
            }
            uTx /= S[j];
            for (label p = 0; p < N; p++)
            {
                y[p] += V(p, j)*uTx;
            }
        }
    }
    x = y;
    forAll(C, c)
    {
        x[nA + c] *= sC[c];
    }

    dLambda = 0;
    forAll(el, p)
    {
        if (!std::isfinite(x[p]))
        {
            return false;
        }
        dLambda[el[p]] = x[p];
    }
    nC.setSize(C.size());
    forAll(C, c)
    {
        if (!std::isfinite(x[nA + c]))
        {
            return false;
        }
        nC[c] = x[nA + c];
    }

    return true;
}


template<class ThermoType>
bool Foam::blastChemicalEquilibrium<ThermoType>::solveTV
(
    const scalar T,
    const scalar V,
    const scalarField& b,
    const boolList& activeEl,
    const boolList& active,
    scalarField& lambda,
    DynamicList<label>& C,
    scalarField& n,
    scalarField& dndT
) const
{
    using constant::thermodynamic::RR;
    using constant::thermodynamic::Pstd;

    const label nE = elements_.size();

    scalarField g(nSpecie_);
    gibbs(T, g);
    const scalar lnc0 = log(V*Pstd/(RR*T));

    scalarField r(nE, 0.0);
    scalarField rc;
    scalarField dLambda(nE, 0.0);
    scalarField nC;
    labelList removed(nSpecie_, 0);

    bool converged = false;
    for (label iter = 0; iter < maxIter_; iter++)
    {
        gasAmounts(lnc0, g, active, lambda, n);

        // Element balance residual and condensed-phase conditions
        r = 0;
        forAll(n, k)
        {
            forAll(r, i)
            {
                r[i] += a_[k][i]*n[k];
            }
        }
        r = b - r;
        rc.setSize(C.size());
        forAll(C, c)
        {
            rc[c] = g[C[c]];
            forAll(lambda, i)
            {
                rc[c] -= a_[C[c]][i]*lambda[i];
            }
        }

        if (!solveKKT(n, b, activeEl, C, r, rc, dLambda, nC))
        {
            if (debug_)
            {
                Info<< "    TV T " << T << " iter " << iter
                    << " singular KKT, C " << C << endl;
            }
            if (C.size())
            {
                C.remove();
                continue;
            }
            return false;
        }

        const scalar dMax = max(mag(dLambda));
        scalar scale = dMax > maxStep_ ? maxStep_/dMax : 1;

        // Backtracking line search on the Lagrangian
        //     L = sum_k n_k(lambda) - b.lambda + sum_c nC_c (a_c.lambda - g_c),
        // convex in lambda, for which the KKT step is a descent direction
        // (dL = -dLambda.H.dLambda)
        {
            auto L = [&](const scalarField& lam, const scalarField& nk)
            {
                scalar l = sum(nk);
                forAll(lam, i)
                {
                    l -= b[i]*lam[i];
                }
                forAll(C, c)
                {
                    scalar x = -g[C[c]];
                    forAll(lam, i)
                    {
                        x += a_[C[c]][i]*lam[i];
                    }
                    l += nC[c]*x;
                }
                return l;
            };

            scalar slope = 0;
            forAll(dLambda, i)
            {
                scalar ri = r[i];
                forAll(C, c)
                {
                    ri -= a_[C[c]][i]*nC[c];
                }
                slope -= ri*dLambda[i];
            }

            const scalar L0 = L(lambda, n);
            scalarField lambdaT(lambda);
            scalarField nT(n);
            for (label ls = 0; ls < 30; ls++)
            {
                lambdaT = lambda + scale*dLambda;
                gasAmounts(lnc0, g, active, lambdaT, nT);
                if
                (
                    slope >= 0
                 || L(lambdaT, nT) <= L0 + 1e-4*scale*slope + 1e-14*mag(L0)
                )
                {
                    break;
                }
                scale *= 0.5;
            }
        }

        lambda += scale*dLambda;

        if (debug_)
        {
            Info<< "    TV T " << T << " iter " << iter << " C " << C
                << " nC " << nC << " dLambda " << dLambda << " lambda "
                << lambda << " r " << r << endl;
        }

        // Remove the condensed species with the most negative amount
        label cMin = -1;
        forAll(nC, c)
        {
            if (nC[c] < 0 && (cMin < 0 || nC[c] < nC[cMin]))
            {
                cMin = c;
            }
        }
        if (cMin >= 0)
        {
            removed[C[cMin]]++;
            for (label c = cMin; c < C.size() - 1; c++)
            {
                C[c] = C[c + 1];
            }
            C.remove();
            continue;
        }

        // Add the condensed species most supersaturated: at any iteration
        // if clearly supersaturated (the gas-phase problem may have no
        // solution, e.g. carbon without oxygen), otherwise near convergence;
        // a species removed before is only added again near convergence
        const bool nearConverged = scale == 1 && dMax < 1e-3;
        label kMax = -1;
        scalar vMax = 0;
        forAll(condensed_, k)
        {
            if (condensed_[k] && active[k] && findIndex(C, k) < 0)
            {
                scalar v = -g[k];
                forAll(lambda, i)
                {
                    v += a_[k][i]*lambda[i];
                }
                const scalar vMin =
                    nearConverged ? 1e-9
                  : removed[k] == 0 ? 1e-3
                  : great;
                if (v > max(vMax, vMin))
                {
                    vMax = v;
                    kMax = k;
                }
            }
        }
        if (kMax >= 0)
        {
            C.append(kMax);
            continue;
        }

        // Converged when the change of every species amount is below the
        // tolerance relative to its amount plus 1e-4 of the total, and the
        // elements balance. Minor species are thus converged to an absolute
        // 1e-14 of the total (with the default tolerance): along directions
        // weakly determined by such species (e.g. the potentials of elements
        // held mainly in a condensed phase) the round-off of the element
        // balance limits the relative accuracy of their amounts
        if (scale == 1)
        {
            const scalar nFloor = 1e-4*(sum(n) + sum(nC));
            scalar dnMax = 0;
            forAll(n, k)
            {
                if (n[k] > 0)
                {
                    scalar x = 0;
                    forAll(dLambda, i)
                    {
                        x += a_[k][i]*dLambda[i];
                    }
                    dnMax = max(dnMax, n[k]*mag(x)/(n[k] + nFloor));
                }
            }
            // Element balance (with the condensed amounts of the KKT
            // solution), relative to the abundance of each element
            scalar rMax = 0;
            forAll(activeEl, i)
            {
                if (activeEl[i])
                {
                    scalar ri = r[i];
                    forAll(C, c)
                    {
                        ri -= a_[C[c]][i]*nC[c];
                    }
                    rMax = max(rMax, mag(ri)/b[i]);
                }
            }

            if (debug_)
            {
                Info<< "    TV dnMax " << dnMax << " rMax " << rMax << endl;
            }
            if (dnMax < tolerance_ && rMax < tolerance_)
            {
                converged = true;
                break;
            }
        }
    }

    if (!converged)
    {
        return false;
    }

    gasAmounts(lnc0, g, active, lambda, n);
    forAll(C, c)
    {
        n[C[c]] = nC[c];
    }

    // Temperature derivative of the composition at fixed V and b:
    // d ln n_k/dT = -1/T - dg_k/dT + a_k.dlambda/dT (gas),
    // a_c.dlambda/dT = dg_c/dT (condensed),
    // sum_k a_k dn_k/dT = 0
    scalarField dgdT(nSpecie_);
    {
        const scalar dT = 1e-5*T;
        scalarField gp(nSpecie_), gm(nSpecie_);
        gibbs(T + dT, gp);
        gibbs(T - dT, gm);
        dgdT = (gp - gm)/(2*dT);
    }
    r = 0;
    forAll(n, k)
    {
        if (!condensed_[k])
        {
            forAll(r, i)
            {
                r[i] += a_[k][i]*n[k]*(1/T + dgdT[k]);
            }
        }
    }
    rc.setSize(C.size());
    forAll(C, c)
    {
        rc[c] = dgdT[C[c]];
    }
    scalarField dndTC;
    if (!solveKKT(n, b, activeEl, C, r, rc, dLambda, dndTC))
    {
        return false;
    }
    forAll(n, k)
    {
        if (!condensed_[k])
        {
            scalar x = -1/T - dgdT[k];
            forAll(dLambda, i)
            {
                x += a_[k][i]*dLambda[i];
            }
            dndT[k] = n[k]*x;
        }
        else
        {
            dndT[k] = 0;
        }
    }
    forAll(C, c)
    {
        dndT[C[c]] = dndTC[c];
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::blastChemicalEquilibrium<ThermoType>::blastChemicalEquilibrium
(
    const speciesTable& species,
    const PtrList<ThermoType>& thermos,
    const dictionary& dict
)
:
    thermos_(thermos),
    nSpecie_(thermos.size()),
    elements_(),
    a_(nSpecie_),
    condensed_(nSpecie_, false),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-10)),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 200)),
    maxStep_(dict.lookupOrDefault<scalar>("maxStep", 2)),
    Tmin_(dict.lookupOrDefault<scalar>("Tmin", 200)),
    Tmax_(dict.lookupOrDefault<scalar>("Tmax", 1e5)),
    debug_(dict.lookupOrDefault<bool>("debug", false))
{
    List<List<Tuple2<word, scalar>>> comp(nSpecie_);
    forAll(species, k)
    {
        const word& name = species[k];
        if (dict.found("elements") && dict.subDict("elements").found(name))
        {
            comp[k] = List<Tuple2<word, scalar>>
            (
                dict.subDict("elements").lookup(name)
            );
        }
        else
        {
            comp[k] = parseFormula(name);
        }
        condensed_[k] = name.find('(') != string::npos;

        forAll(comp[k], j)
        {
            if (findIndex(elements_, comp[k][j].first()) < 0)
            {
                elements_.append(comp[k][j].first());
            }
        }
    }

    if (dict.found("condensedSpecies"))
    {
        const wordList cs(dict.lookup("condensedSpecies"));
        condensed_ = false;
        forAll(cs, i)
        {
            condensed_[species[cs[i]]] = true;
        }
    }

    forAll(a_, k)
    {
        a_[k].setSize(elements_.size(), 0.0);
        forAll(comp[k], j)
        {
            a_[k][findIndex(elements_, comp[k][j].first())] +=
                comp[k][j].second();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::blastChemicalEquilibrium<ThermoType>::elementAbundance
(
    const scalarField& Y
) const
{
    tmp<scalarField> tb(new scalarField(elements_.size(), 0.0));
    scalarField& b = tb.ref();
    forAll(Y, k)
    {
        const scalar nk = max(Y[k], 0)/thermos_[k].W();
        forAll(b, i)
        {
            b[i] += a_[k][i]*nk;
        }
    }
    return tb;
}


template<class ThermoType>
Foam::scalar Foam::blastChemicalEquilibrium<ThermoType>::Ea
(
    const scalarField& Y,
    const scalar rho,
    const scalar T
) const
{
    scalar ea = 0;
    forAll(Y, k)
    {
        ea += Y[k]*(thermos_[k].Es(rho, 0, T) + thermos_[k].Hf());
    }
    return ea;
}


template<class ThermoType>
bool Foam::blastChemicalEquilibrium<ThermoType>::solveTV
(
    const scalar rho,
    const scalar T,
    scalarField& Y
) const
{
    const scalar V = 1/rho;
    const scalarField b(elementAbundance(Y));
    boolList activeEl, active;
    activeSets(b, activeEl, active);

    scalarField n0(nSpecie_);
    forAll(n0, k)
    {
        n0[k] = max(Y[k], 0)/thermos_[k].W();
    }
    scalarField lambda(elements_.size(), 0.0);
    initialLambda(T, V, n0, b, activeEl, active, lambda);

    DynamicList<label> C;
    scalarField n(nSpecie_), dndT(nSpecie_);
    if (!solveTV(T, V, b, activeEl, active, lambda, C, n, dndT))
    {
        return false;
    }

    forAll(Y, k)
    {
        Y[k] = n[k]*thermos_[k].W();
    }
    Y /= sum(Y);
    return true;
}


template<class ThermoType>
bool Foam::blastChemicalEquilibrium<ThermoType>::solveUV
(
    const scalar rho,
    const scalar ea,
    scalarField& Y,
    scalar& T
) const
{
    const scalar V = 1/rho;
    const scalarField b(elementAbundance(Y));
    boolList activeEl, active;
    activeSets(b, activeEl, active);

    scalarField n(nSpecie_), dndT(nSpecie_);
    forAll(n, k)
    {
        n[k] = max(Y[k], 0)/thermos_[k].W();
    }

    T = min(max(T, Tmin_), Tmax_);
    scalarField lambda(elements_.size(), 0.0);
    initialLambda(T, V, n, b, activeEl, active, lambda);
    DynamicList<label> C;

    scalar Tlo = Tmin_, Thi = Tmax_;
    bool bracketLo = false, bracketHi = false;

    // Last temperature with a converged TV solution, and its potentials
    scalar Tgood = -1;
    scalarField lambdaGood(lambda);
    DynamicList<label> Cgood;

    for (label iter = 0; iter < maxIter_; iter++)
    {
        if (!solveTV(T, V, b, activeEl, active, lambda, C, n, dndT))
        {
            // Restart the element potentials from the input composition
            forAll(n, k)
            {
                n[k] = max(Y[k], 0)/thermos_[k].W();
            }
            initialLambda(T, V, n, b, activeEl, active, lambda);
            C.clear();
            if (!solveTV(T, V, b, activeEl, active, lambda, C, n, dndT))
            {
                if (Tgood < 0 || mag(T - Tgood) < 1e-6*Tgood)
                {
                    return false;
                }

                // Retreat towards the last converged temperature
                if (T < Tgood)
                {
                    Tlo = T;
                }
                else
                {
                    Thi = T;
                }
                T = 0.5*(T + Tgood);
                lambda = lambdaGood;
                C = Cgood;
                continue;
            }
        }
        Tgood = T;
        lambdaGood = lambda;
        Cgood = C;

        scalar E = 0, dEdT = 0, cvT = 0;
        forAll(n, k)
        {
            const scalar W = thermos_[k].W();
            const scalar ek = thermos_[k].Es(rho, 0, T) + thermos_[k].Hf();
            const scalar cvk = thermos_[k].Cv(rho, 0, T);
            E += n[k]*W*ek;
            dEdT += n[k]*W*cvk + dndT[k]*W*ek;
            cvT += n[k]*W*cvk*T;
        }
        const scalar f = E - ea;

        if (mag(f) < tolerance_*cvT)
        {
            forAll(Y, k)
            {
                Y[k] = n[k]*thermos_[k].W();
            }
            Y /= sum(Y);
            return true;
        }

        if (f < 0)
        {
            Tlo = T;
            bracketLo = true;
        }
        else
        {
            Thi = T;
            bracketHi = true;
        }

        scalar Tnew = T - f/max(dEdT, small);
        if (dEdT <= 0 || Tnew <= Tlo || Tnew >= Thi)
        {
            if (bracketLo && bracketHi)
            {
                Tnew = 0.5*(Tlo + Thi);
            }
            else if (f < 0)
            {
                Tnew = min(2*T, Thi);
            }
            else
            {
                Tnew = max(0.5*T, Tlo);
            }
        }

        if (mag(Tnew - T) < 1e-14*T)
        {
            break;
        }
        T = Tnew;
    }

    return false;
}


// ************************************************************************* //
