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
#include "varEntry.H"
#include "paramEntry.H"
#include "optEqn.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // argList::addNote
    // (
    //     "Optimises a given case setup by reading scalar entries to modify\n"
    //     "and changing the inputs of the simulation within the given bounds.\n"
    //     "The \"vars\" entry provides the path to the entries to change, and\n"
    //     "the \"errorCriteria\" gives the list of error calculators used to\n"
    //     "determine a cost function.\n"
    //     "\n"
    //     "Entries:\n"
    //     "    Scalar entry: \"path/to/entry\" (lower upper) initial\n"
    //     "    List entry: \"path/to/entry\" (lower upper) initial index\n"
    //     "\n\n"
    //     "Examples:\n\n"
    //     "    gas\n"
    //     "    { \n"
    //     "        ignitor\n"
    //     "        {\n"
    //     "            mass 0.0023;\n"
    //     "        }\n"
    //     "    }\n\n"
    //     "    projectile\n"
    //     "    { \n"
    //     "        resistanceProfile\n"
    //     "        {\n"
    //     "            pressure (1e6 10e7 21e6 35e6);\n"
    //     "        }\n"
    //     "    }\n"
    //     "\n"
    //     "\n"
    //     "    // Method for optimising scalar entries\n"
    //     "    Scalar entry:\n"
    //     "        // Modify the mass entry in the gas/ignitor subDictionary\n"
    //     "        // The lower bound is 1e-5 and the upper bound is 0.001\n"
    //     "        // The initial value is 0.0023 (value above)\n"
    //     "        vars\n"
    //     "        (\n"
    //     "            \"gas/ignitor/mass\" (1e-5 0.001) 0.0023\n"
    //     "        );\n"
    //     "\n"
    //     "\n"
    //     "    // Method for optimising and any entry with multiple scalars\n"
    //     "    // between parenthesis, i.e. ( 0 1 ... N )\n"
    //     "    List entry:\n"
    //     "        // Modify the pressure entry in the projectile/resistanceProfile "
    //     "subDictionary\n"
    //     "        // The lower bound is 20e-6 and the upper bound is 30e6\n"
    //     "        // The initial value is 21e6 (value above)\n"
    //     "        // We are changing index 2 or the 3rd entry "
    //     "(numbering starts from 0)\n"
    //     "        vars\n"
    //     "        (\n"
    //     "            \"projectile/resistanceProfile/pressure\" (20e6 30e6) 21e6 2\n"
    //     "        );\n"
    // );

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
    List<varEntry> variables(optimizationProperties.lookup("variables"));
    if (variables.size() == 1)
    {
        eqnPtr.set(new optEqn1(runTime, optimizationProperties));
    }
    else
    {
        eqnPtr.set(new optEqn(runTime, optimizationProperties));
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

    scalar error = eqn.fX(x, 0);
    //
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
        Info<< results[i] << endl;
    }
    Info<< "Final error = " << error << endl;

    Info<< decrIndent << endl;
    IOobject::writeDivider(Info);

    Info<< nl
        << "Finished" << endl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    return 0;
}


// ************************************************************************* //
