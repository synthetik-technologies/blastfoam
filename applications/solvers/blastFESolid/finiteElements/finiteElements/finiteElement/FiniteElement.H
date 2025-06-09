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

Class
    Foam::FiniteElement

Description
    Templated class for finite element based on element type

SourceFiles
    FiniteElement.C

\*---------------------------------------------------------------------------*/

#ifndef FiniteElement_H
#define FiniteElement_H

#include "finiteElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class FiniteElement Declaration
\*---------------------------------------------------------------------------*/

template<ElementType::Type T>
class FiniteElement
:
    public finiteElement
{
public:

    //- Construct from order and size
    inline FiniteElement(const label order, const label size)
    :
        finiteElement(order, size, T)
    {}

    // Access

        //- Return the number of dimensions
        inline label nDims() const
        {
            return Geometry<T>::nDims;
        }

        //- Return the number of vertices
        inline label nVertices() const
        {
            return Geometry<T>::nVertices;
        }

        // Return number of edges
        inline label nEdges() const
        {
            return Geometry<T>::nEdges;
        }

        //- Return the number of faces
        inline label nFaces() const
        {
            return Geometry<T>::nFaces;
        }

        // Return integration rule for a given order
        virtual const integrationRule& ir
        (
            const label orderi,
            const label orderj,
            const label orderk
        ) const
        {
            return integrationRule::getRule(T, orderi);
        }

};


template<class Element>
class ShellFiniteElement
:
    public Element
{
public:

    //- Runtime type information
    TypeName(Element::typeName_())

    // Constructors
    using Element::Element;

    // Return integration rule for a given order
    virtual const integrationRule& ir
    (
        const label orderi,
        const label orderj,
        const label orderk
    ) const
    {
        return shellIntegrationRule::getRule(this->type_, orderi, orderk);
    }
};


#define addShellFE(Name)                        \
typedef ShellFiniteElement<Name> Name##Shell;   \
defineTemplateTypeNameAndDebug(Name##Shell, 0); \
addNamedTemplateToRunTimeSelectionTable         \
(                                               \
    shellFiniteElement,                         \
    ShellFiniteElement,                         \
    Name,                                       \
    typeOrder,                                  \
    Name                                        \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
