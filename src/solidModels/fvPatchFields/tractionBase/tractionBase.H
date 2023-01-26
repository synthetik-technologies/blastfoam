/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    tractionBase

Description
    Base for all traction boundary conditions

SourceFiles
    tractionBase.C

\*---------------------------------------------------------------------------*/

#ifndef tractionBase_H
#define tractionBase_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class tractionBase Declaration
\*---------------------------------------------------------------------------*/

class tractionBase
{
public:
    //- Static data
    static bool canRelax;

private:
    // Private Data

        const fvPatch& patch_;

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;

        //- Integrated force (neighbour)
        vector force_;

protected:

    //- Update integrated force
    void updateForce();


public:

    //- Runtime type information
    TypeName("solidTraction");


    // Constructors

        //- Construct from patch
        tractionBase(const fvPatch&);

        //- Construct from patch and dictionary
        tractionBase
        (
            const fvPatch&,
            const dictionary&,
            const bool mustRead = true
        );

        //- Construct by copy
        tractionBase
        (
            const tractionBase&,
            const fvPatch&
        );

        //- Construct by mapping given tractionBase onto a new patch
        tractionBase
        (
            const tractionBase&,
            const fvPatch&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tractionBase(const tractionBase&) = delete;


    //- Destructor
    virtual ~tractionBase()
    {}



    // Member functions

        // Access

            //- Return integrated force
            inline const vector& force() const
            {
                return force_;
            }

            //- Return traction
            inline const vectorField& traction() const
            {
                return traction_;
            }

            //- Access traction
            inline vectorField& traction()
            {
                return traction_;
            }

            //- Return pressure
            inline const scalarField& pressure() const
            {
                return pressure_;
            }

            //- Access pressure
            inline scalarField& pressure()
            {
                return pressure_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const tractionBase&,
                const labelList&
            );

        //- Update pressure and traction
        virtual bool updateFields() = 0;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
