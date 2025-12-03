/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

Class
    Foam::masterSystem

Description
    Class to handle a group of systems that need to be updated once per group
SourceFiles
    masterSystem.C

\*---------------------------------------------------------------------------*/

#ifndef masterSystemList_H
#define masterSystemList_H

#include "regIOobject.H"
#include "UPtrList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class masterSystem;

/*---------------------------------------------------------------------------*\
                 Class masterSystemList Declaration
\*---------------------------------------------------------------------------*/

class masterSystemList
:
    public regIOobject,
    public UPtrList<masterSystem>
{
public:

    TypeName("masterSystemList");

    // Constructor
    masterSystemList(const objectRegistry& obr);

    //- Lookup or construct masterSystemList
    static masterSystemList& New(const objectRegistry& obr);


    //- Destructor
    ~masterSystemList();

    //- Disallow default bitwise copy construct
    masterSystemList(const masterSystem&) = delete;


    // Member Functions

        //- Add phase to list of kinetic theory phases.
        virtual void addSystem(masterSystem&);

        //- Initialize the model after everything is constructed
        virtual void initialize();

        //- Correct total volume fraction
        virtual void correctAlpha();

        //- Update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

       //- Explicit post update
        virtual void solveExplicit();

        //- Save explicit changes
        virtual void storeExplicit();

        //- Implicit post update
        virtual void solveImplicit();

        //- Post update
        virtual void postUpdate();

        //- Remove stored fields
        virtual void clear();

        virtual bool writeData(Ostream&) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const masterSystemList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
