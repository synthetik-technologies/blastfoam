/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::searchableAnnulusCylinder

Description
    Searching on annulus cylinder

SourceFiles
    searchableAnnulusCylinder.C

\*---------------------------------------------------------------------------*/

#ifndef searchableAnnulusCylinder_H
#define searchableAnnulusCylinder_H

#include "treeBoundBox.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                Class searchableAnnulusCylinder Declaration
\*---------------------------------------------------------------------------*/

class searchableAnnulusCylinder
:
    public searchableSurface
{
private:
    // Private Member Data

        //- 'left' point
        const point point1_;

        //- 'right' point
        const point point2_;

        //- Length of vector point2-point1
        const scalar magDir_;

        //- Normalised vector point2-point1
        const vector unitDir_;

        //- Inner radius
        const scalar innerRadius_;

        //- Outer radius
        const scalar outerRadius_;

        //- Names of regions
        mutable wordList regions_;


    // Private Member Functions

        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        //- Find nearest point on cylinder.
        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        scalar radius2(const point& pt) const;

        //- Find intersection with cylinder
        label findLineAll
        (
            const point& start,
            const point& end,
            pointIndexHit& nearO,
            pointIndexHit& nearI,
            pointIndexHit& farI,
            pointIndexHit& farO
        ) const;

        //- Return the boundBox of the cylinder
        boundBox calcBounds() const;


public:

    //- Runtime type information
    TypeName("searchableAnnulusCylinder");


    // Constructors

        //- Construct from components
        searchableAnnulusCylinder
        (
            const IOobject& io,
            const point&,
            const point&,
            const scalar innerRadius,
            const scalar outerRadius
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableAnnulusCylinder
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        searchableAnnulusCylinder(const searchableAnnulusCylinder&) = delete;


    //- Destructor
    virtual ~searchableAnnulusCylinder();


    // Member Functions

        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return true;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 1;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const
        {
            NotImplemented;

            return false;
        }


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
                return false;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchableAnnulusCylinder&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
