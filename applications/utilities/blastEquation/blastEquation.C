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

#include "Time.H"



namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


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


enum EqnType
{
    SINGLE,
    UNI,
    MULTI,
    UNKNOWN
};


EqnType getEqnType(const string& fxStr)
{
    label i = fxStr.find("x[");
    if (i < 0)
    {
        return SINGLE;
    }
    return UNI;
}


EqnType getEqnType(const argList& args)
{
    if (args.optionFound("P"))
    {
        return SINGLE;
    }
    return getEqnType(args.optionRead<string>("fx"));
}


EqnType getEqnType(const dictionary& dict)
{
    if (dict.found("P"))
    {
        return SINGLE;
    }
    return getEqnType(dict.lookup<verbatimString>("fx_code"));
}


label calcNVar(const string& fxStr)
{
    IStringStream is(fxStr);
    label maxI = -1;
    while (is.good())
    {
        char t(readChar(is));
        if (t == token::BEGIN_SQR)
        {
            label i(readLabel(is));
            maxI = max(i, maxI);
        }
    }
    return maxI + 1;
}

word readXs(Istream& is)
{
    OStringStream os;
    word str;
    while (is.good())
    {
        token t(is);
        if (t.good())
        {
            os << t << token::SPACE;
        }
    }
    return os.str();
}


List<scalar> readSingleX(Istream& is)
{
    DynamicList<scalar> xs;
    token t(is);
    while (is.good())
    {
        is >> t;
        if (t.isNumber())
        {
            xs.append(t.number());
        }
    }
    return move(xs);
}


List<scalarList> readMultiX(Istream& is)
{
    DynamicList<scalarList> xs;
    DynamicList<scalar> x;
    token t(is); // read first (
    bool inList = false;

    is >> t;
    // Check for multiple inputs
    if (t.isPunctuation())
    {
        if (t.pToken() == token::BEGIN_LIST)
        {
            inList = true;
        }
    }
    else
    {
        inList = true;
        is.putBack(t);
    }
    while (is.good())
    {
        is >> t;

        if (t.isPunctuation())
        {
            if (t.pToken() == token::BEGIN_LIST)
            {
                inList = true;
            }
            else if (t.pToken() == token::END_LIST && inList)
            {
                inList = false;
                xs.append(x);
                x.clear();
            }
        }
        else if (t.isNumber())
        {
            if (inList)
            {
                x.append(t.number());
            }
        }
    }
    return move(xs);
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
        if (args.optionFound("dfdx"))
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
        dict.set("P", P);
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
    else
    {
        FatalErrorInFunction
            << "For scalar equations, either fx of P (polynomial coefficients)"
            << " must be provided" << endl
            << abort(FatalError);
    }
    dict.set("name", args.optionLookupOrDefault<word>("name", "f"));

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
    }
    else if (nDerivatives == 0)
    {
        nDerivatives++;
    }

    if (!dict.found("d2fdx2_code"))
    {
        dict.set("d2fdx2_code", string("return 0.0;"));
    }
    else if (nDerivatives == 1)
    {
        nDerivatives++;
    }

    if (!dict.found("d3fdx3_code"))
    {
        dict.set("d3fdx3_code", string("return 0.0;"));
    }
    else if (nDerivatives == 2)
    {
        nDerivatives++;
    }

    return nDerivatives;
}


void setEquationSolverDict
(
    const argList& args,
    const dictionary& funcDict,
    const label nDerivatives,
    dictionary& dict
)
{
    if (args.optionFound("eval"))
    {
        dictionary evaluationDict;
        evaluationDict.set("x", readXs((args.optionLookup("eval")())));
        dict.set("evaluationCoeffs", evaluationDict);
        dict.set("evaluate", true);
    }

    if (args.optionFound("findRoots") || args.optionFound("findAllRoots"))
    {
        dictionary rootsDict;
        rootsDict.set
        (
            "solver",
            args.optionLookupOrDefault<word>
            (
                "rootSolver",
                (nDerivatives > 0 ? "NewtonRaphson" : "secant")
            )
        );
        if (args.optionFound("findAllRoots"))
        {
            rootsDict.set("findAllRoots", true);
            Pair<scalar> rootBounds
            (
                args.optionRead<Pair<scalar>>("findAllRoots")
            );
            rootsDict.set("bounds", rootBounds);
        }
        else
        {
            Pair<scalar> rootBounds
            (
                args.optionRead<Pair<scalar>>("findRoots")
            );
            rootsDict.set("bounds", rootBounds);
            rootsDict.set
            (
                "x0",
                0.5*(rootBounds.first() + rootBounds.second())
            );
        }
        dict.set("rootCoeffs", rootsDict);
        dict.set("findRoots", true);
    }

    if (args.optionFound("findPolyRoots") && funcDict.found("P"))
    {
        dict.set("findRoots", true);
    }

    if (args.optionFound("integrate"))
    {
        dictionary integrationDict;
        integrationDict.set
        (
            "bounds",
            args.optionRead<Pair<scalar>>("integrate")
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

        List<scalar> xs(args.optionRead<List<scalar>>("minimize"));
        minimizationDict.set
        (
            "bounds",
            Pair<scalar>(xs[0], xs[1])
        );
        if (xs.size() > 2)
        {
            minimizationDict.set("x0", xs[2]);
        }

        minimizationDict.set
        (
            "solver",
            args.optionLookupOrDefault<word>
            (
                "minimizer",
                "bisection"
            )
        );

        dict.set("minimizationCoeffs", minimizationDict);
        dict.set("minimize", true);
    }
}


string splitEquations(const word& name, Istream& is)
{
    DynamicList<string> strs;
    autoPtr<OStringStream> os;
    token t;

    while (is.good())
    {
        if (!os.valid())
        {
            os.set(new OStringStream());
            os() << name << "[" << strs.size() << "] = ";
        }
        is >> t;
        if (t.good())
        {
            if (t.isPunctuation())
            {
                if (t.pToken() == token::END_STATEMENT)
                {
                    os() << t;
                    strs.append(os->str());
                    os.clear();
                    continue;
                }
            }
            os() << t << token::SPACE;
        }
    }
    if (os.valid())
    {
        if (os().str().size())
        {
            os() << token::END_STATEMENT;
            strs.append(os->str());
        }
    }

    string str;
    forAll(strs, i)
    {
        str += strs[i];
    }
    return str;

}
void setUnivariateFuncDict
(
    const argList& args,
    dictionary& dict
)
{
    dict.set("fx_code", "return " + args.optionRead<string>("fx") +";");
    dict.set("eqnString", args.optionRead<string>("fx"));
    if (args.optionFound("dfdx"))
    {
        dict.set
        (
            "dfdx_code",
            splitEquations("dfdx", (args.optionLookup("dfdx"))())
        );
    }

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
        nVar = calcNVar(args.optionRead<string>("fx"));
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
    if (!dict.found("dfdx_code"))
    {
        dict.set
        (
            "dfdx_code",
            string("UnivariateEquation<scalar>::dfdX(x, li, dfdx);")
        );
    }
    return 1;
}


void setUnivariateSolverDict
(
    const argList& args,
    const dictionary& funcDict,
    const label nDerivatives,
    dictionary& dict
)
{
    if (args.optionFound("eval"))
    {
        dictionary evaluationDict;
        evaluationDict.set("x", readXs((args.optionLookup("eval"))()));
        dict.set("evaluationCoeffs", evaluationDict);
        dict.set("evaluate", true);
    }

    if (args.optionFound("integrate"))
    {
        dictionary integrationDict;
        integrationDict.set
        (
            "bounds",
            args.optionRead<Pair<scalarList>>("integrate")
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
        integrationDict.set
        (
            "nNodes",
            labelList(funcDict.lookup<label>("nVar"), 3)
        );
        dict.set("integrate", true);

        dict.set("integrationCoeffs", integrationDict);
    }

    if (args.optionFound("minimize"))
    {
        dictionary minimizationDict;
        List<scalarList> xs(args.optionRead<List<scalarList>>("minimize"));
        minimizationDict.set
        (
            "bounds",
            Pair<scalarList>(xs[0], xs[1])
        );
        if (xs.size() > 2)
        {
            minimizationDict.set("x0", xs[2]);
        }

        minimizationDict.set
        (
            "solver",
            args.optionLookupOrDefault<word>
            (
                "minimizer",
                "gradientDescent"
            )
        );
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
    argList::addOption("dict", "Operations dictionary");

    argList::addBoolOption("clean", "Remove dynamicCode library");
    argList::addBoolOption("noClean", "No not remove dynamicCode library");
    argList::addOption("name", "Function name");

    argList::addOption("fx", "Function");
    argList::addOption("dfdx", "Derivative function");
    argList::addOption("d2fdx2", "Second derivative function");
    argList::addOption("d3fdx3", "Third derivative function");
    argList::addOption("P", "Polynomial coeffs");

    argList::addOption("eval", "Evaluate the function at the given values");

    argList::addOption("findRoots", "Find nearest root");
    argList::addOption("findAllRoots", "Find all roots");
    argList::addOption("rootSolver", "Root solver type");
    argList::addBoolOption("findPolyRoots", "Find roots of a polynomial");

    argList::addOption("minimize", "Find function minimum");
    argList::addOption("minimizer", "Minimization scheme");

    argList::addOption("integrate", "Integrate function");
    argList::addOption("integrator", "Integration scheme");

    argList::addBoolOption("time", "Construct time");

    #include "setRootCase.H"

    autoPtr<Time> runTimePtr;
    if (args.optionFound("time"))
    {
        Foam::Info<< "Create time\n" << Foam::endl;

        runTimePtr.set(new Time(Foam::Time::controlDictName, args));
    }
    else
    {
        fileName path(cwd());
        runTimePtr.set
        (
            new Time
            (
                fileName(path.path()),
                fileName(path.name())
            )
        );
    }
    Time& runTime = runTimePtr();


    autoPtr<dictionary> funcDictPtr;
    autoPtr<dictionary> dictPtr;

    label nDerivatives = -1;
    EqnType eqnType = UNKNOWN;
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
        eqnType = getEqnType(funcDictPtr());

        if (funcDictPtr->found("P"))
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
    }
    else
    {
        eqnType = getEqnType(args);

        clean = !args.optionFound("noClean") || clean;
        funcDictPtr.set(new dictionary());
        if (eqnType == UNI)
        {
            setUnivariateFuncDict(args, funcDictPtr());
        }
        else
        {
            setEquationFuncDict(args, funcDictPtr(), P);
        }
    }
    if (eqnType == SINGLE)
    {
        nDerivatives = nEquationDerivatives(funcDictPtr());
    }
    else if (eqnType == UNI)
    {
        nDerivatives = nUnivariateEquationDerivatives(funcDictPtr());
    }

    if (!funcDictPtr->found("nDerivatives"))
    {
        funcDictPtr->set("nDerivatives", nDerivatives);
    }

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
        if (eqnType == UNI)
        {
            setUnivariateSolverDict
            (
                args,
                funcDictPtr(),
                nDerivatives,
                dictPtr()
            );
        }
        else if (eqnType == SINGLE)
        {
            setEquationSolverDict
            (
                args,
                funcDictPtr(),
                nDerivatives,
                dictPtr()
            );
        }

    }
    const dictionary& dict = dictPtr();


    if (eqnType == UNI)
    {
        const word name = funcDictPtr->lookupOrDefault<word>("name", "f");

        CodedUnivariateEquation<scalar> eqn(runTime, funcDictPtr());
        Info << nl << endl;

        if (eqn.eqnString() != "undefined")
        {
            Info<< "************************************" << nl
                << name << "(x) = " << eqn.eqnString() << nl
                << "************************************" << nl
                << endl;
        }

        if (dict.lookupOrDefault("evaluate", false))
        {
            const dictionary& evaluationDict(dict.subDict("evaluationCoeffs"));

            Info<< "************************************" << nl
                << "Evaluating function" << nl
                << "************************************" << endl;

            List<scalarList> xs
            (
                readMultiX(evaluationDict.lookup("x"))
            );
            forAll(xs, i)
            {
                Info<<"f(" << print(xs[i]) << ") = "
                    << eqn.fX(xs[i], 0) << endl;
                if (nDerivatives > 0)
                {
                    scalarList dfdx(xs[i].size());
                    eqn.dfdX(xs[i], 0, dfdx);
                    Info<<"dfdx(" << print(xs[i]) << ") = " << dfdx << endl;
                }
                Info<< endl;
            }
        }

        if (dict.lookupOrDefault("integrate", false))
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
            Info<<"integral from " << bounds[0] << " to " << bounds[1] << " = "
                << integrator->integrate(bounds[0], bounds[1], 0) << nl
                << nl << endl;
        }

        if (dict.lookupOrDefault("minimize", false))
        {
            const dictionary& minimizationDict
            (
                dict.subDict("minimizationCoeffs")
            );

            Info<< "************************************" << nl
                << "Minimizing function" << nl
                << "************************************" << endl;

            autoPtr<minimizationScheme> minimizer
            (
                minimizationScheme::New(eqn, minimizationDict)
            );
            const scalarList lower(eqn.lowerLimits());
            const scalarList upper(eqn.upperLimits());
            if (minimizationDict.found("bounds"))
            {
                Pair<scalarList> limits(minimizationDict.lookup("bounds"));
                eqn.setLowerLimits(limits.first());
                eqn.setUpperLimits(limits.second());
            }

            scalarList minimum
            (
                minimizationDict.found("x0")
              ? minimizer->solve
                (
                    minimizationDict.lookup<scalarList>("x0"),
                    0
                )
              : minimizer->solve()
            );
            Info<< "Minimum = " << minimum << nl << endl;

            eqn.setLowerLimits(lower);
            eqn.setUpperLimits(upper);
        }
    }
    else if (eqnType == SINGLE)
    {
        const word name = funcDictPtr->lookupOrDefault<word>("name", "f");

        CodedEquation<scalar> eqn(runTime, funcDictPtr());
        Info << nl << endl;

        if (eqn.name() != "undefined")
        {
            Info<< "************************************" << nl
                << name << "(x) = " << eqn.name() << nl
                << "************************************" << nl
                << endl;
        }

        if (dict.lookupOrDefault("evaluate", false))
        {
            const dictionary& evaluationDict(dict.subDict("evaluationCoeffs"));

            Info<< "************************************" << nl
                << "Evaluating function" << nl
                << "************************************" << endl;

            List<scalar> xs(readSingleX(evaluationDict.lookup("x")));
            forAll(xs, i)
            {
                Info<<"f(" << xs[i] << ") = " << eqn.fx(xs[i], 0) << endl;
                if (nDerivatives > 0)
                {
                    Info<<"dfdx(" << xs[i] << ") = "
                        << eqn.dfdx(xs[i], 0) << endl;
                }
                if (nDerivatives > 1)
                {
                    Info<<"d2fdx2(" << xs[i] << ") = "
                        << eqn.d2fdx2(xs[i], 0) << endl;
                }
                if (nDerivatives > 2)
                {
                    Info<<"d3fdx3(" << xs[i] << ") = "
                        << eqn.d3fdx3(xs[i], 0) << endl;
                }
            }

            Info<< endl;
        }

        if (dict.lookupOrDefault("findRoots", false))
        {
            Info<< "************************************" << nl
                << "Finding function roots" << nl
                << "************************************" << endl;
            if (P.size())
            {
                polynomialRoots PRoots(P);
                Info<< "Real roots = " << PRoots.rootsRe() << nl
                    << "Imaginary roots = " << PRoots.rootsIm() << nl
                    << endl;
            }
            else
            {
                const dictionary& rootsDict(dict.subDict("rootCoeffs"));

                autoPtr<univariateRootSolver> rootSolver
                (
                    univariateRootSolver::New(eqn, rootsDict)
                );
                Pair<scalar> bounds(rootsDict.lookup("bounds"));

                if (rootsDict.lookupOrDefault("findAllRoots", false))
                {
                    label nIntervals
                    (
                        rootsDict.lookupOrDefault("nIntervals", 123)
                    );
                    scalarList roots
                    (
                        rootSolver->solveAll
                        (
                            bounds[0],
                            bounds[1],
                            0,
                            nIntervals
                        )
                    );
                    Info<< "Roots = " << roots << nl << endl;
                }
                else
                {
                    scalar x(rootsDict.lookup<scalar>("x0"));
                    scalar root
                    (
                        rootSolver->solve(x, bounds[0], bounds[1], 0.0)
                    );

                    Info<< "Root = " << root << nl << endl;
                }
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
            Info<<"integral from " << bounds[0] << " to " << bounds[1] << " = "
                << integrator->integrate(bounds[0], bounds[1], 0) << nl
                << endl;
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

            autoPtr<minimizationScheme> minimizer
            (
                minimizationScheme::New(eqn, minimizationDict)
            );
            const scalar lower(eqn.lower());
            const scalar upper(eqn.upper());
            if (minimizationDict.found("bounds"))
            {
                Pair<scalar> limits(minimizationDict.lookup("bounds"));
                eqn.setLower(limits.first());
                eqn.setUpper(limits.second());
            }

            scalar minimum
            (
                minimizationDict.found("x0")
              ? minimizer->solve
                (
                    {minimizationDict.lookup<scalar>("x0")},
                    0
                )()[0]
              : minimizer->solve()()[0]
            );
            Info<< "Minimum = " << minimum << nl << endl;

            eqn.setLower(lower);
            eqn.setUpper(upper);
        }
    }

    if (clean)
    {
        rmDir("./dynamicCode");
    }

    Info<< "done" << endl;
    return 0;
}


// ************************************************************************* //
