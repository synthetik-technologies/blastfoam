/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::equationBase

Description
    Base class for all equation types that holds the equation as an optional
    string

\*---------------------------------------------------------------------------*/

#ifndef equationBase_H
#define equationBase_H

#include "dictionary.H"
#include "objectRegistry.H"
#include "OFstream.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class equationBase Declaration
\*---------------------------------------------------------------------------*/

class equationBase
{
    // Private data

        //- Optional pointer to object registry
        UautoPtr<const objectRegistry> obrPtr_;

        //- Optional log file
        bool log_;
        fileName logFile_;
        bool append_;
        mutable autoPtr<OFstream> logPtr_;


protected:
    // Protected data

        //- Equation as a string
        word eqnString_;

        //- Merge list of strings
        string mergeStrings(const List<string>& eqns) const;

        //- Return the object registry
        const objectRegistry& obr() const;


public:

    // Constructors

        //- Construct null
        equationBase();

        //- Construct from equation string
        equationBase(const string& eqnString);

        //- Construct from multiple equation strings
        equationBase(const List<string>& eqnStrings);

        //- Construct from equation string and dictionary
        equationBase(const string& eqnString, const dictionary& dict);

        //- Construct from multiple equation strings and dictionary
        equationBase(const List<string>& eqnStrings, const dictionary& dict);

        //- Construct from a dictionary
        equationBase(const dictionary& dict);


    //- Destructor
    virtual ~equationBase();


    // Member functions

        //- Return the equation as a string
        virtual const word& eqnString() const
        {
            return eqnString_;
        }

        //- Access the equation string
        virtual word& eqnString()
        {
            return eqnString_;
        }

        //- Set the object registry
        virtual void setObr(const objectRegistry& obr)
        {
            obrPtr_.reset(&obr);
        }

        //- Return if logging is wanted
        bool log() const
        {
            return log_;
        }

        //- Set the log name
        void setLog(const fileName&, const bool log = true);

        //- Set the log name
        void setLog(const bool log);

        //- Set the append flag for logging
        void setAppend(const bool append = true)
        {
            append_ = append;
        }

        //- Return the number of derivatives
        virtual label nDerivatives() const = 0;

        //- Return the log stream
        OFstream& logStream() const;

        //- Read logging information
        virtual void read(const dictionary&);
};

#define defineEquationTypeNameAndDebug(Eqn, Debug, Type) \
defineTemplateTypeNameAndDebug(Eqn<Type>, Debug);

#define defineEquationTypeName(Eqn, Type) \
defineTemplateTypeName(Eqn<Type>);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
