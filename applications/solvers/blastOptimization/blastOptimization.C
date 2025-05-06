/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022 Synthetik Applied Technologies
     \\/     M anipulation  |
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

Application
    blastOptimisation

Description
    Optimization solver for varying IB inputs to restrict specified field values


\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "IOdictionary.H"
#include "minimizationScheme.H"
#include "variableEntry.H"
#include "parameterEntry.H"
#include "optimizationEquation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Optimises a given case setup by reading scalar entries to modify\n"
        "and changing the inputs of the simulation within the given bounds.\n"
        "The \"variables\" entry provides the path to the entries to change, and\n"
        "the \"errors\" gives the list of error calculators used to\n"
        "determine a cost function. The \"optError\" functionObject should be\n"
        "used in the actual simulation so the error is written\n"
        "\n"
        "Entries:\n"
        "    Scalar entry: \"path/to/entry\" (lower upper) initial\n"
        "    List entry: \"path/to/entry\" [index] (lower upper) initial \n"
        "\n\n"
        "Examples:\n\n"
        "    dict00\n"
        "    { \n"
        "        dict01\n"
        "        {\n"
        "            val 0.0023;\n"
        "        }\n"
        "    }\n\n"
        "    dict10\n"
        "    { \n"
        "        dict11\n"
        "        {\n"
        "            lst (0 1 2);\n"
        "        }\n"
        "    }\n"
        "\n"
        "\n"
        "    // Method for optimising scalar entries\n"
        "    Scalar entry:\n"
        "        // Modify the mass entry in the gas/ignitor subDictionary\n"
        "        // The lower bound is 0 and the upper bound is 1\n"
        "        // The initial value is 0.5 (value above)\n"
        "        variables\n"
        "        (\n"
        "            \"dict00/dict01/val\" (0 1) 0.5\n"
        "        );\n"
        "\n"
        "\n"
        "    // Method for optimising and any entry with multiple scalars\n"
        "    // between parenthesis, i.e. ( 0 1 ... N )\n"
        "    List entry:\n"
        "        // Modify the entry in the dict10/dict11 subDictionary\n"
        "        // The lower bound is 0 and the upper bound is 10\n"
        "        // The initial value is 5 (value above)\n"
        "        // We are changing index 2 or the 3rd entry (numbering starts from 0)\n"
        "        variable\n"
        "        (\n"
        "            \"dict10/dict11/lst\" [2] (0 10) 5\n"
        "        );\n"
    );
    argList::addBoolOption("restart", "Used cached results to restart simulation");

    #include "setRootCase.H"
    #include "createTime.H"

    //- Creating an IOdictionary
    IOdictionary optimizationProperties
    (
        IOobject
        (
            "optimizationProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Set up optimization
    autoPtr<scalarUnivariateEquation> eqnPtr;
    List<variableEntry> variables(optimizationProperties.lookup("variables"));
    if (variables.size() == 1)
    {
        eqnPtr.set
        (
            new optimizationEquation1
            (
                runTime,
                optimizationProperties,
                args.optionFound("restart")
            )
        );
    }
    else
    {
        eqnPtr.set
        (
            new optimizationEquation
            (
                runTime,
                optimizationProperties,
                args.optionFound("restart")
            )
        );
    }
    scalarUnivariateEquation& eqn = eqnPtr();

    IOobject::writeDivider(Info);
    Info<< "Optimization variables:" << incrIndent << endl;
    List<scalar> x(variables.size());
    forAll(variables, i)
    {
        Info<< variables[i] << endl;
        x[i] = variables[i].value();
    }
    decrIndent(Info);
    IOobject::writeDivider(Info);

    // Set debug level so each step is printed
    minimizationScheme::debug = 3;
    autoPtr<minimizationScheme> solverPtr
    (
        minimizationScheme::New(eqn, optimizationProperties)
    );

    Info<< endl;
    IOobject::writeDivider(Info);
    Info<< endl;

    List<scalar> results = solverPtr->solve(x);

    // Run optimized case
    scalar error = eqn.fX(results, 0);

    label nSteps = max(solverPtr->nSteps(), 1);

    Info<< nl;
    IOobject::writeDivider(Info);
    Info<< nl;

    if (solverPtr->converged())
    {
        Info<< "Converged in " << nSteps << " iterations" << endl;
    }
    else
    {
        Warning
            << "**** Did not converged in " << nSteps << " iterations" << endl;
    }
    Info<< "Optimized values are: " << incrIndent << endl;
    forAll(variables, i)
    {
        Info<< indent
            << fileName(variables[i].file())
            << "|" << fileName(variables[i].path())
            << " = " << results[i] << endl;
    }
    Info<< nl << decrIndent
        << "Final error = " << error << nl << endl;


    IOobject::writeDivider(Info);

    Info<< nl
        << "Finished" << endl
        << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    return 0;
}


// ************************************************************************* //
