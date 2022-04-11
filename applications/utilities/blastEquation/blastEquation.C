/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-10-2019  Jeff Heylmun:   Added refinement to setFields utility
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

Description
    Set values on a selected set of cells/patchfaces through a dictionary and
    refines using hexRef method.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"

#include "EquationsFwd.H"
#include "univariateRootSolver.H"
#include "IntegratorsFwd.H"
#include "univariateMinimizationScheme.H"
#include "polynomialRoots.H"

#include "UnivariateEquationsFwd.H"
#include "rootSolver.H"
#include "MultivariateIntegratorsFwd.H"
#include "minimizationScheme.H"

#include "CodedEquation.H"
#include "CodedUnivariateEquation.H"
#include "CodedMultivariateEquation.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef CodedEquation<scalar> scalarCodedEquation;
defineTemplateTypeNameAndDebug(scalarCodedEquation, 0);
addTemplateToRunTimeSelectionTable
(
    scalarEquation,
    CodedEquation,
    scalar,
    dictionary
);

typedef CodedUnivariateEquation<scalar> scalarCodedUnivariateEquation;
defineTemplateTypeNameAndDebug(scalarCodedUnivariateEquation, 0);
addTemplateToRunTimeSelectionTable
(
    scalarUnivariateEquation,
    CodedUnivariateEquation,
    scalar,
    dictionary
);

typedef CodedMultivariateEquation<scalar> scalarCodedMultivariateEquation;
defineTemplateTypeNameAndDebug(scalarCodedMultivariateEquation, 0);
addTemplateToRunTimeSelectionTable
(
    scalarMultivariateEquation,
    CodedMultivariateEquation,
    scalar,
    dictionary
);


word print(const scalarList& x)
{
    OStringStream os;
    for (label i = 0; i < x.size() - 1; i++)
    {
        os << x[i] << ", ";
    }
    os << x.last();
    return os.str();
}


void setEquationFuncDict
(
    const argList& args,
    dictionary& dict,
    scalarList& P
)
{
    if (args.optionFound("fx"))
    {
        dict.set("eqnString", args.optionRead<string>("fx"));
        dict.set("fx_code", "return " + args.optionRead<string>("fx") +";");
        if (!args.optionFound("dfdx"))
        {
            dict.set
            (
                "dfdx_code",
                "return " + args.optionRead<string>("dfdx") + ";"
            );
        }
        if (args.optionFound("d2fdx2"))
        {
            dict.set
            (
                "d2fdx2_code",
                "return " + args.optionRead<string>("d2fdx2") + ";"
            );
        }
        if (args.optionFound("d3fdx3"))
        {
            dict.set
            (
                "d3fdx3_code",
                "return " + args.optionRead<string>("d3fdx3") + ";"
            );
        }
    }
    else if (args.optionFound("P"))
    {
        P = args.optionRead<List<scalar>>("P");
        dict.set("fx_code", "return " + polynomialRoots::polyName(P) + ";");
        dict.set("eqnString", string(polynomialRoots::polyName(P)));

        List<scalar> dfdx(polynomialRoots::derivative(P));
        dict.set
        (
            "dfdx_code",
            "return " + polynomialRoots::polyName(dfdx) + ";"
        );

        List<scalar> d2fdx2(polynomialRoots::derivative(dfdx));
        dict.set
        (
            "d2fdx2_code",
            "return " + polynomialRoots::polyName(d2fdx2) + ";"
        );

        List<scalar> d3fdx3(polynomialRoots::derivative(d2fdx2));
        dict.set
        (
            "d3fdx3_code",
            "return " + polynomialRoots::polyName(d3fdx3) + ";"
        );
    }
    word name
    (
        string(args.optionLookupOrDefault<word>("name", "f"))
    );
    dict.set("name", string(name));

    Pair<scalar> bounds(-great, great);
    if (args.optionFound("bounds"))
    {
        bounds = args.optionRead<Pair<scalar>>("bounds");
    }
    dict.set("bounds", bounds);
    dict.set("lowerBound", bounds[0]);
    dict.set("upperBound", bounds[1]);
}


label nEquationDerivatives(dictionary& dict)
{
    label nDerivatives = 0;
    if (!dict.found("dfdx_code"))
    {
        dict.set("dfdx_code", string("return 0.0;"));
        if (nDerivatives == 0)
        {
            nDerivatives++;
        }
    }
    if (!dict.found("d2fdx2_code"))
    {
        dict.set("d2fdx2_code", string("return 0.0;"));
        if (nDerivatives == 1)
        {
            nDerivatives++;
        }
    }
    if (!dict.found("d3fdx3_code"))
    {
        dict.set("d3fdx3_code", string("return 0.0;"));
        if (nDerivatives == 2)
        {
            nDerivatives++;
        }
    }
    return nDerivatives;
}


void setEquationSolverDict
(
    const argList& args,
    const label nDerivatives,
    dictionary& dict
)
{
    if (args.optionFound("eval"))
    {
        dictionary evaluationDict;
        if (args.optionFound("x"))
        {
            evaluationDict.set("x", args.optionRead<scalar>("x"));
        }
        if (args.optionFound("xs"))
        {
            evaluationDict.set("xs", args.optionRead<scalarList>("xs"));
        }

        dict.set("evaluationCoeffs", evaluationDict);
        dict.set("evaluate", true);
    }

    if (args.optionFound("root"))
    {
        dictionary rootDict;
        Pair<scalar> rootBounds
        (
            args.optionLookupOrDefault<Pair<scalar>>
            (
                "rootBounds",
                dict.lookup("bounds")
            )
        );
        rootDict.set("bounds", rootBounds);
        rootDict.set
        (
            "solver",
            args.optionLookupOrDefault<word>
            (
                "rootSolver",
                (nDerivatives > 0 ? "NewtonRaphson" : "secant")
            )
        );
        rootDict.set
        (
            "x",
            args.optionLookupOrDefault<scalar>
            (
                "x",
                0.5*(rootBounds.first() + rootBounds.second())
            )
        );

        dict.set("rootCoeffs", rootDict);
        dict.set("findRoot", true);
    }
    if (args.optionFound("allRoots"))
    {
        dictionary rootsDict;
        rootsDict.set
        (
            "rootBounds",
            args.optionRead<Pair<scalar>>("rootBounds")
        );
        rootsDict.set
        (
            "solver",
            args.optionLookupOrDefault<word>
            (
                "rootSolver",
                nDerivatives > 0 ? "NewtonRaphson" : "secant"
            )
        );

        dict.set("rootsCoeffs", rootsDict);
        dict.set("findAllRoots", true);
    }

    if (args.optionFound("integrate"))
    {
        dictionary integrationDict;
        integrationDict.set
        (
            "bounds",
            args.optionRead<Pair<scalar>>("integrationBounds")
        );
        integrationDict.set
        (
            "integrator",
            args.optionLookupOrDefault<word>
            (
                "integrator",
                "Simpson13"
            )
        );

        dict.set("integrationCoeffs", integrationDict);
        dict.set("integrate", true);

    }

    if (args.optionFound("minimize"))
    {
        dictionary minimizationDict;

        minimizationDict.set
        (
            "solver",
            args.optionLookupOrDefault<word>
            (
                "minimizer",
                "bisection"
            )
        );
        dict.set("x0", args.optionRead<scalar>("x"));

        dict.set("minimizationCoeffs", minimizationDict);
        dict.set("minimize", true);
    }
}


void setUnivariateFuncDict
(
    const argList& args,
    dictionary& dict
)
{
    dict.set("fx_code", "return " + args.optionRead<string>("fx") +";");
    dict.set("eqnString", args.optionRead<string>("fx"));

    dict.set("name", string(args.optionLookupOrDefault<word>("name", "f")));

    Pair<scalarList> bounds;
    label nVar = -1;
    if (args.optionFound("bounds"))
    {
        bounds = args.optionRead<Pair<scalarList>>("bounds");
        nVar = bounds[0].size();
    }
    else
    {
        nVar = args.optionRead<label>("nVar");
        bounds[0].setSize(nVar, -great);
        bounds[1].setSize(nVar, great);
    }
    dict.set("bounds", bounds);
    dict.set("lowerBounds", bounds[0]);
    dict.set("upperBounds", bounds[1]);
    dict.set("nVar", nVar);
}


label nUnivariateEquationDerivatives(dictionary& dict)
{
    label nDerivatives = 0;
    if (!dict.found("dfdx_code"))
    {
        dict.set
        (
            "dfdx_code",
            string("UnivariateEquation<scalar>::dfdX(x, li, dfdx);")
        );
        if (nDerivatives == 0)
        {
            nDerivatives++;
        }
    }
    if (!dict.found("d2fd2_code"))
    {
        dict.set
        (
            "d2fdx2_code",
            string("UnivariateEquation<scalar>::d2fdX2(x, li, d2fdx2);")
        );
        if (nDerivatives == 1)
        {
            nDerivatives++;
        }
    }
    if (!dict.found("d3fdx3_code"))
    {
        dict.set
        (
            "d3fdx3_code",
            string("UnivariateEquation<scalar>::d3fdX3(x, li, d3fdx3);")
        );
        if (nDerivatives == 2)
        {
            nDerivatives++;
        }
    }
    return nDerivatives;
}


void setUnivariateSolverDict
(
    const argList& args,
    const label nVar,
    const label nDerivatives,
    dictionary& dict
)
{
    if (args.optionFound("eval"))
    {
        dictionary evaluationDict;
        if (args.optionFound("x"))
        {
            evaluationDict.set("x", args.optionRead<scalarList>("x"));
        }
        if (args.optionFound("xs"))
        {
            evaluationDict.set("xs", args.optionRead<List<scalarList>>("xs"));
        }

        dict.set("evaluationCoeffs", evaluationDict);
        dict.set("evaluate", true);
    }

    if (args.optionFound("integrate"))
    {
        dictionary integrationDict;
        integrationDict.set
        (
            "bounds",
            args.optionRead<Pair<scalarList>>("integrationBounds")
        );
        integrationDict.set
        (
            "integrator",
            args.optionLookupOrDefault<word>
            (
                "integrator",
                "Gaussian"
            )
        );
        integrationDict.set("nNodes", labelList(nVar, 3));
        dict.set("integrate", true);

        dict.set("integrationCoeffs", integrationDict);
    }

    if (args.optionFound("minimize"))
    {
        dictionary minimizationDict;

        minimizationDict.set
        (
            "solver",
            args.optionLookupOrDefault<word>
            (
                "minimizer",
                "NelderMead"
            )
        );
        dict.set("x0", args.optionRead<scalarList>("x"));

        dict.set("minimizationCoeffs", minimizationDict);
        dict.set("minimize", true);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;
int main(int argc, char *argv[])
{
    argList::addOption("func", "Function dictionary");
    argList::addBoolOption("clean", "Remove dynamicCode library");
    argList::addBoolOption("noClean", "No not remove dynamicCode library");
    argList::addOption("name", "Function name");

    argList::addBoolOption("multi", "Use multiple inputs for function");
    argList::addOption("nVar", "numberOfVariables");

    argList::addOption("x", "Evaluation point");
    argList::addOption("fx", "Function");
    argList::addOption("dfdx", "Derivative function");
    argList::addOption("d2fdx2", "Second derivative function");
    argList::addOption("d3fdx3", "Third derivative function");
    argList::addOption("P", "Polynomial coeffs");
    argList::addOption("bounds", "bounds");

    argList::addBoolOption("eval", "Evaluate the function");
    argList::addOption("x", "Value to evaluate at");
    argList::addOption("xs", "Values to evaluate at");

    argList::addBoolOption("root", "Find nearest root");
    argList::addBoolOption("allRoots", "Find all roots");
    argList::addOption("rootSolver", "Root solver type");
    argList::addOption("rootBounds", "Bounds of the roots");

    argList::addBoolOption("minimize", "Find function minimum");
    argList::addOption("minimizer", "Minimization scheme");

    argList::addBoolOption("integrate", "Integrate function");
    argList::addOption("integrator", "Integration scheme");
    argList::addOption("integrationBounds", "Bounds of integration");

    #include "addDictOption.H"

    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    autoPtr<dictionary> funcDictPtr;
    autoPtr<dictionary> dictPtr;

    label nDerivatives = -1;
    bool multi = false;
    List<scalar> P;

    bool clean = args.optionFound("clean");

    if (args.optionFound("func"))
    {
        funcDictPtr.set
        (
            new dictionary
            (
                (IFstream(args.optionRead<fileName>("func")))()
            )
        );
        multi =
            args.optionFound("multi")
          || funcDictPtr->lookupOrDefault("multi", false);

        if (funcDictPtr->found("P") && !multi)
        {
            P = funcDictPtr->lookup<scalarList>("P");
            funcDictPtr->set
            (
                "fx_code",
                "return " + polynomialRoots::polyName(P) + ";"
            );
            funcDictPtr->set("eqnString", string(polynomialRoots::polyName(P)));

            List<scalar> dfdx(polynomialRoots::derivative(P));
            funcDictPtr->set
            (
                "dfdx_code",
                "return " + polynomialRoots::polyName(dfdx) + ";"
            );

            List<scalar> d2fdx2(polynomialRoots::derivative(dfdx));
            funcDictPtr->set
            (
                "d2fdx2_code",
                "return " + polynomialRoots::polyName(d2fdx2) + ";"
            );

            List<scalar> d3fdx3(polynomialRoots::derivative(d2fdx2));
            funcDictPtr->set
            (
                "d3fdx3_code",
                "return " + polynomialRoots::polyName(d3fdx3) + ";"
            );
        }
        nDerivatives =
            multi
          ? nUnivariateEquationDerivatives(funcDictPtr())
          : nEquationDerivatives(funcDictPtr());

        if (!funcDictPtr->found("nDerivatives"))
        {
            funcDictPtr->set("nDerivatives", nDerivatives);
        }
    }
    else
    {
        clean = !args.optionFound("noClean") || clean;
        funcDictPtr.set(new dictionary());
        multi = args.optionFound("multi");
        if (multi)
        {
            setUnivariateFuncDict(args, funcDictPtr());
            nDerivatives = nUnivariateEquationDerivatives(funcDictPtr());
        }
        else
        {
            setEquationFuncDict(args, funcDictPtr(), P);
            nDerivatives = nEquationDerivatives(funcDictPtr());
        }
    }
    funcDictPtr->set("nDerivatives", nDerivatives);

    if (args.optionFound("dict"))
    {
        dictPtr.set
        (
            new dictionary
            (
                (IFstream(args.optionRead<fileName>("dict")))()
            )
        );
    }
    else
    {
        dictPtr.set(new dictionary());
        if (multi)
        {
            setUnivariateSolverDict
            (
                args,
                funcDictPtr->lookup<label>("nVar"),
                nDerivatives,
                dictPtr()
            );
        }
        else
        {
            setEquationSolverDict(args, nDerivatives, dictPtr());
        }

    }
    const dictionary& dict = dictPtr();


    if (multi)
    {
        const word name = dict.lookupOrDefault<word>("name", "f");

        CodedUnivariateEquation<scalar> eqn(funcDictPtr());
        Info << nl << endl;

        if (eqn.name() != "undefined")
        {
            Info<< "************************************" << nl
                << name << "(x)=" << eqn.name() << nl
                << "************************************" << nl
                << endl;
        }

        if (dict.isDict("evaluationCoeffs"))
        {
            const dictionary& evaluationDict(dict.subDict("evaluationCoeffs"));

            Info<< "************************************" << nl
                << "Evaluating function" << nl
                << "************************************" << endl;

            bool evaluated = false;
            if (evaluationDict.found("x"))
            {
                scalarList x(evaluationDict.lookup("x"));
                scalarList dfdx(x.size());
                Info<<"f(" << print(x) << ") = " << eqn.fX(x, 0) << endl;
                eqn.dfdX(x, 0, dfdx);
                evaluated = true;
            }
            if (dict.found("xs"))
            {
                List<scalarList> xs(evaluationDict.lookup("xs"));
                forAll(xs, i)
                {
                    Info<<"f(" << xs[i] << ") = " << eqn.fX(xs[i], 0) << endl;
                }
                evaluated = true;
            }
            if (!evaluated)
            {
                WarningInFunction
                    << "Please provide \"x\" or \"xs\" to evaluate a function"
                    << endl;
            }
            Info<< endl;
        }

        if (dict.isDict("integrationCoeffs"))
        {
            const dictionary& integrationDict(dict.subDict("integrationCoeffs"));

            Info<< "************************************" << nl
                << "Integrating function" << nl
                << "************************************" << endl;


            Pair<scalarList> bounds(integrationDict.lookup("bounds"));
            autoPtr<MultivariateIntegrator<scalar>> integrator
            (
                MultivariateIntegrator<scalar>::New(eqn, integrationDict)
            );
            Info<<"integral_{" << bounds[0] << "}^{" << bounds[1] << "}: "
                << integrator->integrate(bounds[0], bounds[1], 0) << nl
                << integrator->nIntervals() << " intervals"
                << nl << endl;
        }

        if (dict.isDict("minimizationCoeffs"))
        {
            const dictionary& minimizationDict
            (
                dict.subDict("minimizationCoeffs")
            );

            Info<< "************************************" << nl
                << "Minimizing function" << nl
                << "************************************" << endl;

            minimizationScheme::debug = 2;
            autoPtr<minimizationScheme> minimizer
            (
                minimizationScheme::New(eqn, minimizationDict)
            );
            minimizer->solve
            (
                minimizationDict.lookup<scalarList>("x0"),
                0
            );
            Info<< endl;
        }
    }
    else
    {
        const word name = dict.lookupOrDefault<word>("name", "f");

        CodedEquation<scalar> eqn(funcDictPtr());
        Info << nl << endl;

        if (eqn.name() != "undefined")
        {
            Info<< "************************************" << nl
                << name << "(x)=" << eqn.name() << nl
                << "************************************" << nl
                << endl;
        }

        if (dict.isDict("evaluationCoeffs"))
        {
            const dictionary& evaluationDict(dict.subDict("evaluationCoeffs"));

            Info<< "************************************" << nl
                << "Evaluating function" << nl
                << "************************************" << endl;

            bool evaluated = false;
            if (evaluationDict.found("x"))
            {
                scalar x = evaluationDict.lookup<scalar>("x");
                Info<<"f(" << x << ") = " << eqn.fx(x, 0) << endl;
                evaluated = true;
            }
            if (dict.found("xs"))
            {
                List<scalar> xs(evaluationDict.lookup<List<scalar>>("xs"));
                forAll(xs, i)
                {
                    Info<<"f(" << xs[i] << ") = " << eqn.fx(xs[i], 0) << endl;
                }
                evaluated = true;
            }
            if (!evaluated)
            {
                WarningInFunction
                    << "Please provide \"x\" or \"xs\" to evaluate a function"
                    << endl;
            }
            Info<< endl;
        }

        if (dict.isDict("rootCoeffs"))
        {
            const dictionary& rootDict(dict.subDict("rootCoeffs"));

            Info<< "************************************" << nl
                << "Finding function root" << nl
                << "************************************" << endl;

            autoPtr<univariateRootSolver> rootSolver
            (
                univariateRootSolver::New(eqn, rootDict)
            );
            Pair<scalar> bounds(rootDict.lookup("bounds"));
            scalar x(rootDict.lookup<scalar>("x0"));

            Info<<"Root: "
                << rootSolver->solve(x, bounds[0], bounds[1], 0.0)
                << nl << endl;
        }

        if (dict.isDict("rootsCoeffs"))
        {
            const dictionary& rootsDict(dict.subDict("rootsCoeffs"));

            Info<< "************************************" << nl
                << "Finding function roots" << nl
                << "************************************" << endl;

            if (P.size() > 1)
            {
                polynomialRoots PRoots(P);
                Info<< "Real roots: " << PRoots.rootsRe() << nl
                    << "Imaginary roots: " << PRoots.rootsIm() << nl
                    << endl;
            }
            else
            {
                autoPtr<univariateRootSolver> rootSolver
                (
                    univariateRootSolver::New(eqn, rootsDict)
                );

                Pair<scalar> bounds(rootsDict.lookup("bounds"));
                label nIntervals
                (
                    rootsDict.lookupOrDefault("nIntervals", 123)
                );
                Info<< "Root: "
                    <<  rootSolver->solveAll
                        (
                            bounds[0],
                            bounds[1],
                            0,
                            nIntervals
                        )
                    << nl << endl;
            }
        }

        if (dict.isDict("integrationCoeffs"))
        {
            const dictionary& integrationDict(dict.subDict("integrationCoeffs"));

            Info<< "************************************" << nl
                << "Integrating function" << nl
                << "************************************" << endl;


            Pair<scalar> bounds(integrationDict.lookup("bounds"));
            autoPtr<Integrator<scalar>> integrator
            (
                Integrator<scalar>::New(eqn, dictPtr())
            );
            Info<<"integral_{" << bounds[0] << "}^{" << bounds[1] << "}: "
                << integrator->integrate(bounds[0], bounds[1], 0) << nl
                << integrator->nIntervals() << " intervals"
                << nl << endl;
        }

        if (dict.isDict("minimizationCoeffs"))
        {
            const dictionary& minimizationDict
            (
                dict.subDict("minimizationCoeffs")
            );

            Info<< "************************************" << nl
                << "Minimizing function" << nl
                << "************************************" << endl;

            minimizationScheme::debug = 2;
            autoPtr<univariateMinimizationScheme> minimizer
            (
                univariateMinimizationScheme::New(eqn, minimizationDict)
            );

            minimizer->solve
            (
                minimizationDict.lookup<scalar>("x0"),
                0
            );
            Info<< endl;
        }
    }

    if (args.optionFound("clean"))
    {
        rmDir("./dynamicCode");
    }

    Info<< "done" << endl;
    return 0;
}


// ************************************************************************* //
