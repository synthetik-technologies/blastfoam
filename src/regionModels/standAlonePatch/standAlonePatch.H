/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Typedef
    Foam::standAlonePatch

Description
    standAlonePatch is a primitive patch that holds its point and face list.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef standAlonePatch_H
#define standAlonePatch_H

#include "PrimitivePatch.H"
#include "face.H"
#include "List.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class standAlonePatch
:
    public PrimitivePatch<faceList, pointField>
{
    //- Non-const reference to the points
    pointField& pointsRef_;


public:
        //- Construct from components
        standAlonePatch
        (
            const faceList& faces,
            const pointField& points
        );

        //- Move constructor from components
        standAlonePatch
        (
            faceList&& faces,
            pointField&& points
        );

        //- Move constructor from components
        standAlonePatch
        (
            faceList&& faces,
            List<point>&& points
        );

        //- Construct from components, reuse storage
        standAlonePatch
        (
            faceList& faces,
            pointField& points,
            const bool reuse
        );

        //- Copy constructor
        standAlonePatch(const standAlonePatch&);

        //- Move constructor
        standAlonePatch(standAlonePatch&&);

        //- Construct and return a clone
        virtual autoPtr<PrimitivePatch<faceList, pointField>> clone() const
        {
            return autoPtr<PrimitivePatch<faceList, pointField>>
            (
                new standAlonePatch(*this)
            );
        }

        //- Construct from Istream
        standAlonePatch(Istream&, const pointField&);


    //- Destructor
    virtual ~standAlonePatch();


    // Member Functions

        using PrimitivePatch<faceList, pointField>::points;


        // Edit

            //- Correct patch after moving points
            virtual void movePoints(const pointField&);


    // Member Operators

        //- Assignment operator
        void operator=(const standAlonePatch&);

        //- Assignment operator
        template<class PrimitivePatchRef>
        void operator=(const PrimitivePatchRef& pp)
        {
            clearOut();
            faceList::operator=(pp);
            pointsRef_ = pp.points();
        }

        //- Move assignment operator
        void operator=(standAlonePatch&&);

        //- Assignment operator
        template<class PrimitivePatchRef>
        void operator=(PrimitivePatchRef&& pp)
        {
            clearOut();
            faceList::operator=(move(pp));
            pointsRef_ = move(pp.points());
        }

    void writeVTK(const fileName&) const;

    standAlonePatch createGlobalPatch() const;

    template<class Patch>
    static standAlonePatch createGlobalPatch(const Patch& patch);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "standAlonePatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
