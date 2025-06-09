/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
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

\*---------------------------------------------------------------------------*/

#include "finiteElement.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(finiteElement, 0);
    defineRunTimeSelectionTable(finiteElement, type);
    defineRunTimeSelectionTable(finiteElement, typeOrder);
    defineRunTimeSelectionMap(finiteElement, msh);

    defineTypeNameAndDebug(shellFiniteElement, 0);
    defineRunTimeSelectionTable(shellFiniteElement, type);
    defineRunTimeSelectionTable(shellFiniteElement, typeOrder);
//     defineRunTimeSelectionMap(shellFiniteElement, msh);
}

Foam::finiteElement::finiteElementTable
    Foam::finiteElement::finiteElements;
Foam::List<Foam::List<Foam::label>> Foam::finiteElement::Binomials;

Foam::shellFiniteElement::finiteElementTable
    Foam::shellFiniteElement::finiteElements;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::List<Foam::label>& Foam::finiteElement::Binomial(const label p)
{
    if (Binomials.size() <= p)
    {
        Binomials.setSize(p+1);
    }
    if (!Binomials[p].size())
    {
        List<label>& binomial = Binomials[p];
        binomial.setSize(p+1, 0);
        forAll(binomial, i)
        {
            binomial[0] = 1;
            label binij = 1;
            label binij1 = 1;
            for (label j = 1; j < i; j++)
            {
                binij = binomial[j];
                binomial[j] = binomial[j] + binij1;
                binij1 = binij;
            }
            binomial[i] = 1;
        }
    }
    return Binomials[p];
}


void Foam::finiteElement::calcBinomialTerms
(
    const label p,
    const scalar x,
    const scalar y,
    List<label>& u
)
{
    u.setSize(p+1);
    if (p == 0)
    {
        return;
    }

    label i = 0;
    const labelList b(Binomial(p));
    scalar z = x;

    for (i = 1; i < p; i++)
    {
        u[i] = b[i]*z;
        z *= x;
    }

    u[p] = z;
    z = y;
    for  (i--; i > 0; i--)
    {
        u[i] *= z;
        z *= y;
    }
    u[0] = z;
}

void Foam::finiteElement::calcBinomialTerms
(
    const label p,
    const scalar x,
    const scalar y,
    List<label>& u,
    List<label>& du
)
{
    u.setSize(p+1);
    du.setSize(p+1);
    if (p == 0)
    {
        return;
    }

    const scalar xpy = x + y;
    const scalar ptx = scalar(p)*x;
    label i = 0;
    const labelList b(Binomial(p));
    scalar z = 1;

    for (i = 1; i < p; i++)
    {
        du[i] = b[i]*z*(scalar(i)*xpy - ptx);
        z *= x;
        u[i] = b[i]*z;
    }
    du[p] = p*z;
    u[p] = z;
    z = 1;
    for  (i--; i > 0; i--)
    {
        du[i] *= z;
        z *= y;
        u[i] *= z;
    }
    du[0] = -p*z;
    u[0] = z*y;
}


Foam::labelList Foam::finiteElement::dofMap(const label nDim, const label o)
{
    const label o1 = o+1;
    label I = 0;

    List<label> map(pow(o1, nDim));

    if (nDim == 1)
    {
        map[0] = I++;
        map[o] = I++;

        for (label i = 1; i < o; i++) map[i] = I++;

    }
    else if (nDim == 2)
    {
        auto index = [o1](const label i, const label j)
        {
            return i + o1*j;
        };
        map[index(0, 0)] = I++;
        map[index(o, 0)] = I++;
        map[index(o, o)] = I++;
        map[index(0, o)] = I++;

        // (0, 1)
        for (label i = 1; i < o; i++) map[index(i, 0)] = I++;

        // (1, 2)
        for (label j = 1; j < o; j++) map[index(o, j)] = I++;

        // (2, 3)
        for (label i = 1; i < o; i++) map[index(o-i, o)] = I++;

        // (3, 0)
        for (label j = 1; j < o; j++) map[index(0, o-j)] = I++;

        // (0, 1, 2, 3)
        for (label j = 1; j < o; j++)
            for (label i = 1; i < o; i++)
                map[index(i, j)] = I++;

    }
    else if (nDim == 3)
    {
        auto index = [o1](const label i, const label j, const label k)
        {
            return i + o1*(j + o1*k);
        };

        map[index(0, 0, 0)] = I++;
        map[index(o, 0, 0)] = I++;
        map[index(o, o, 0)] = I++;
        map[index(0, o, 0)] = I++;
        map[index(0, 0, o)] = I++;
        map[index(o, 0, o)] = I++;
        map[index(o, o, o)] = I++;
        map[index(0, o, o)] = I++;

        // (0, 1)
        for (label i = 1; i < o; i++) map[index(i, 0, 0)] = I++;

        // (1, 2)
        for (label j = 1; j < o; j++) map[index(o, j, 0)] = I++;

        // (3, 2)
        for (label i = 1; i < o; i++) map[index(i, o, 0)] = I++;

        // (0, 3)
        for (label j = 1; j < o; j++) map[index(0, j, 0)] = I++;

        // (4, 5)
        for (label i = 1; i < o; i++) map[index(i, 0, o)] = I++;

        // (5, 6)
        for (label j = 1; j < o; j++) map[index(o, j, o)] = I++;

        // (7, 6)
        for (label i = 1; i < o; i++) map[index(i, o, o)] = I++;

        // (4, 7)
        for (label j = 1; j < o; j++) map[index(0, j, o)] = I++;

        // (0, 4)
        for (label k = 1; k < o; k++) map[index(0, 0, k)] = I++;

        // (1, 5)
        for (label k = 1; k < o; k++) map[index(o, 0, k)] = I++;

        // (2, 6)
        for (label k = 1; k < o; k++) map[index(o, o, k)] = I++;

        // (3, 7)
        for (label k = 1; k < o; k++) map[index(0, o, k)] = I++;

        // (3, 2, 1, 0)
        for (label j = 1; j < o; j++)
            for (label i = 1; i < o; i++)
                map[index(i, o-j, 0)] = I++;

        // (0, 1, 5, 4)
        for (label i = 1; i < o; i++)
            for (label k = 1; k < o; k++)
                map[index(i, 0, k)] = I++;

        // (1, 2, 6, 5)
        for (label k = 1; k < o; k++)
            for (label j = 1; j < o; j++)
                map[index(o, j, k)] = I++;

        // (2, 3, 7, 6)
        for (label k = 1; k < o; k++)
            for (label i = 1; i < o; i++)
                map[index(o-i, o, k)] = I++;

        // (3, 0, 4, 7)
        for (label k = 1; k < o; k++)
            for (label j = 1; j < o; j++)
                map[index(0, o-j, k)] = I++;

        // (4, 5, 6, 7)
        for (label j = 1; j < o; j++)
            for (label i = 1; i < o; i++)
                map[index(i, j, o)] = I++;


        // internal
        for (label k = 1; k < o; k++)
            for (label j = 1; j < o; j++)
                for (label i = 1; i < o; i++)
                    map[index(i, j, k)] = I++;
    }
    return map;
}


const Foam::finiteElement* Foam::finiteElement::getRefFiniteElement
(
    const ElementType::Type et,
    const label order
)
{
    const word type(ElementType::cellModelNames[et]);
    Tuple2<word, label> key(type, order);
    if (!finiteElements.found(key))
    {
        typename typeConstructorTable::iterator cstrIter =
            typeConstructorTablePtr_->find(type);

        if (cstrIter == typeConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown finiteElement type " << type << nl
                << "Valid finiteElement types are:" << nl
                << typeConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        autoPtr<finiteElement> fePtr(cstrIter()(order));
        finiteElements.insert(key, fePtr.ptr());
    }
    return finiteElements[key];
}


const Foam::finiteElement* Foam::finiteElement::getRefFiniteElement
(
    const word& type,
    const label order
)
{
    Tuple2<word, label> key(type, order);
    if (!finiteElements.found(key))
    {
        typename typeConstructorTable::iterator cstrIter =
            typeConstructorTablePtr_->find(type);

        if (cstrIter == typeConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown finiteElement type " << type << nl
                << "Valid finiteElement types are:" << nl
                << typeConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        autoPtr<finiteElement> fePtr(cstrIter()(order));
        finiteElements.insert(key, fePtr.ptr());
    }
    return finiteElements[key];
}

const Foam::finiteElement* Foam::shellFiniteElement::getRefFiniteElement
(
    const ElementType::Type et,
    const label order
)
{
    const word type(ElementType::cellModelNames[et]);
    Tuple2<word, label> key(type, order);
    if (!finiteElements.found(key))
    {
        typename typeConstructorTable::iterator cstrIter =
            typeConstructorTablePtr_->find(type);

        if (cstrIter == typeConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown finiteElement type " << type << nl
                << "Valid finiteElement types are:" << nl
                << typeConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        autoPtr<finiteElement> fePtr(cstrIter()(order));
        finiteElements.insert(key, fePtr.ptr());
    }
    return finiteElements[key];
}


const Foam::finiteElement* Foam::shellFiniteElement::getRefFiniteElement
(
    const word& type,
    const label order
)
{
    Tuple2<word, label> key(type, order);
    if (!finiteElements.found(key))
    {
        typename typeConstructorTable::iterator cstrIter =
            typeConstructorTablePtr_->find(type);

        if (cstrIter == typeConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown finiteElement type " << type << nl
                << "Valid finiteElement types are:" << nl
                << typeConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        autoPtr<finiteElement> fePtr(cstrIter()(order));
        finiteElements.insert(key, fePtr.ptr());
    }
    return finiteElements[key];
}

namespace Foam
{
template<template<class> class ListT>
labelList _matchNodes
(
    const ListT<vector>& localNodes,
    const List<vector>& globalNodes,
    const labelList& indices,
    const scalar tol
)
{
    labelList mappedIndices(localNodes.size(), -1);
    const UIndirectList<vector> nodes(globalNodes, indices);

    if (localNodes.size() < 10)
    {
        forAll(localNodes, i)
        {
            label ni = -1;
            const vector& n = localNodes[i];
            forAll(nodes, j)
            {
                if (magSqr(n - nodes[j]) < tol)
                {
                    ni = j;
                    break;
                }
            }

            if (ni == -1)
            {
                FatalErrorInFunction
                    << "Unmatched node. Hanging nodes are not supported"
                    << endl
                    << abort(FatalError);
            }

            mappedIndices[i] = indices[ni];

        }
    }
    else
    {
        labelHashSet nodesToAdd(identityMap(indices.size()));
        forAll(localNodes, i)
        {
            const vector& n = localNodes[i];
            label ni = -1;
            forAllConstIter(labelHashSet, nodesToAdd, iter)
            {
                if (magSqr(n - nodes[iter.key()]) < tol)
                {
                    ni = iter.key();
                    break;
                }
            }

            if (ni == -1)
            {
                FatalErrorInFunction
                    << "Unmatched node. Hanging nodes are not supported"
                    << endl
                    << abort(FatalError);
            }

            mappedIndices[i] = indices[ni];
            nodesToAdd.unset(ni);
        }
    }

    return mappedIndices;
}
}

Foam::labelList Foam::finiteElement::matchNodes
(
    const UList<vector>& localNodes,
    const List<vector>& globalNodes,
    const labelList& indices,
    const scalar tol
)
{
    return _matchNodes<UList>(localNodes, globalNodes, indices, tol);
}


Foam::labelList Foam::finiteElement::matchNodes
(
    const UIndirectList<vector>& localNodes,
    const List<vector>& globalNodes,
    const labelList& indices,
    const scalar tol
)
{
    return _matchNodes<UIndirectList>(localNodes, globalNodes, indices, tol);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElement::finiteElement
(
    const label order,
    const label size,
    const ElementType::Type type
)
:
    order_(order),
    ir_(integrationRule::getRule(type, 2*order)),
    nodes_(size, Zero),
    type_(type)
{}

Foam::shellFiniteElement::shellFiniteElement()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElement::~finiteElement()
{}


Foam::shellFiniteElement::~shellFiniteElement()
{}

// ************************************************************************* //

