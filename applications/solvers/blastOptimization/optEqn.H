/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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

Class
    Foam::levelSetModel

Description
    Equation for determining error from a given solver

SourceFiles
    optEqn.C

\*---------------------------------------------------------------------------*/
#ifndef optEqn_H
#define optEqn_H

#include "UnivariateEquationsFwd.H"
#include "EquationsFwd.H"
#include "varEntry.H"
#include "paramEntry.H"
#include "commandEntry.H"

namespace Foam
{

class optEqnBase
{
protected:

    // List of entry tokens for variables
    mutable List<varEntry> variables_;

    // List of different configurations
    mutable List<List<paramEntry>> configurations_;

    //- Commands to run
    List<commandEntry> commands_;

    //- Output file
    fileName results_;

    //- Saved list of variables and result error
    mutable label savedIndex_;
    List<List<scalar>> savedVars_;
    List<scalar> savedErrors_;

    //- Read log file
    void readLogFile(const fileName& file, const label nVar);

    //- Return log file name for a given command number and configuration
    string logFile(const label cmd, const label config) const;

    //- Read limits from variable entries
    List<scalar> limits(const label cmpt) const;

    scalar run() const;


public:
    //Constructor
    optEqnBase(const dictionary& dict);

    //Destructor
    ~optEqnBase();

        //- Return variables
        inline const List<varEntry>& vars() const
        {
            return variables_;
        }
};

class optEqn1
:
    public optEqnBase,
    public ScalarEquation
{
public:
    //Constructor
    optEqn1(const Time& runTime, const dictionary& dict, const bool restart);

    //Destructor
    ~optEqn1();


        //Scalar Function
        virtual scalar fx(const scalar x, const label li) const;

        //- Return the number of derivatives
        virtual label nDerivatives() const
        {
            return 0;
        }
};

class optEqn
:
    public optEqnBase,
    public ScalarUnivariateEquation
{
public:
    //Constructor
    optEqn(const Time& runTime, const dictionary& dict, const bool restart);

    //Destructor
    ~optEqn();

        //Scalar Function
        virtual scalar fX(const VarType& x, const label li) const;

        //- Return the number of derivatives
        virtual label nDerivatives() const
        {
            return 0;
        }
};
}



#endif
