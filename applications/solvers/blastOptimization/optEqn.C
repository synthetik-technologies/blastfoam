#include "optEqn.H"
#include "Time.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "IOmanip.H"

Foam::List<Foam::scalar>
Foam::optEqnBase::limits(const label cmpt) const
{
    List<scalar> lim(variables_.size());
    forAll(variables_, i)
    {
        lim[i] = variables_[i].bounds()[cmpt];
    }
    return lim;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optEqnBase::optEqnBase(const dictionary& dict)
:
    variables_(dict.lookup("variables")),
    configurations_
    (
        dict.lookupOrDefault("configurations", List<List<paramEntry>>(1))
    ),
    commands_(dict.lookup("commands")),
    results_(dict.lookup("results"))
{
    IOobject::writeDivider(Info);
    Info<< "Commands " << incrIndent << endl;
    forAll(commands_, i)
    {
        Info<< indent << commands_[i] << endl;
    }
    Info<< decrIndent << endl;

    if (configurations_.size() && configurations_[0].size())
    {
        IOobject::writeDivider(Info);
        Info<< "Optimization configurations:" << incrIndent << endl;
        forAll(configurations_, i)
        {
            Info<< indent << "configuration " << i << incrIndent << endl;
            forAll(configurations_[i], j)
            {
                Info<< indent << configurations_[i][j] << endl;
            }
            decrIndent(Info);
        }
        Info<< decrIndent << endl;
    }
    if (!configurations_.size())
    {
        configurations_.setSize(1);
    }
}


Foam::optEqn1::optEqn1(const Time& runTime, const dictionary& dict)
:
    optEqnBase(dict),
    ScalarEquation
    (
        variables_[0].bounds()[0],
        variables_[0].bounds()[1],
        dict
    )
{
    ScalarEquation::setLog
    (
        runTime.globalPath() / runTime.globalCaseName() + ".evals",
        this->log()
    );
    ScalarEquation::read(dict);
}


Foam::optEqn::optEqn(const Time& runTime, const dictionary& dict)
:
    optEqnBase(dict),
    ScalarUnivariateEquation
    (
        limits(0),
        limits(1),
        dict
    )
{
    ScalarUnivariateEquation::setLog
    (
        runTime.globalPath() / runTime.globalCaseName() + ".evals",
        this->log()
    );
    ScalarUnivariateEquation::read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optEqnBase::~optEqnBase()
{}


Foam::optEqn1::~optEqn1()
{}


Foam::optEqn::~optEqn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::string Foam::optEqnBase::logFile
(
    const label cmd,
    const label config
) const
{
    if (configurations_.size() > 1)
    {
        return
            "log."
          + commands_[cmd].command()
          + "." + Foam::name(config);
    }
    return "log." + commands_[cmd].command();
}


Foam::scalar Foam::optEqnBase::run() const
{
    scalar totalError = 0.0;

    // Loop all parameters and run commands
    forAll(configurations_, configi)
    {
        // Info<< "Running configuration " << configi << ": " << endl;
        string setEnvCmd
        (
            "export BLAST_CONFIG_NO=" + Foam::name(configi)
          + " && export N_BLAST_CONFIGS=" + Foam::name(configurations_.size())
          + " && "
        );
        incrIndent(Info);

        forAll(configurations_[configi], i)
        {
            configurations_[configi][i].set();
        }

        forAllConstIter
        (
            HashTable<dictionary>,
            optEntry::dictionaries,
            iter
        )
        {
            OFstream os(iter.key());
            iter().write(os, false);
        }

        forAll(commands_, i)
        {
            fileName file(logFile(i, configi));

            // Run command, setting configuration variables
            int status = system((setEnvCmd + commands_[i] + " > " + file).c_str());

            // Check return status
            if (status == 2)
            {
                std::exit(2);
            }
            else if (status)
            {
                FatalError
                    << "[FAILED] command: "<< commands_[i] << endl
                    << "check " << file << " for more details" << endl
                    << exit(FatalError);
            }
        }

        // Read in results
        IFstream is(results_);
        if (!is.good())
        {
            FatalErrorInFunction
                << results_ << " does not exist" << endl
                << abort(FatalError);
        }


        HashTable<Pair<scalar>> errors;
        word errorName;
        scalar errorValue, value;
        while (is.good())
        {
            is >> errorName >> value >> errorValue;
            errors.insert(errorName, {value, errorValue});
            // if (errorName != "totalError")
            // {
            //     Info<< errorName << ": "
            //         << "value = " << value
            //         << ", error = " << errorValue << endl;
            // }
        }

        // Info<< "Total error = " << errors["totalError"][1] << endl;
        totalError += errors["totalError"][1];
    }

    return sqrt(totalError);
}


Foam::scalar Foam::optEqn1::fx(const scalar x, const label li) const
{
    // Set to initial values
    variables_[0].set(x);

    scalar error = run();

    if (this->log())
    {
        const unsigned int w = IOstream::defaultPrecision() + 7;
        this->logStream() << setf(ios_base::left);

        static bool writeHeader = true;
        if (writeHeader)
        {
            writeHeader = false;
            this->logStream()
                << setw(w) << "# x" << token::SPACE
                << setw(w) << "error" << endl;
        }
        this->logStream()
            << setw(w) << x << ' ' << setw(w) << error << endl;
    }
    return error;
}


Foam::scalar Foam::optEqn::fX(const VarType& x, const label li) const
{
    // Set current variables
    forAll(variables_, i)
    {
        //Set to initial values
        variables_[i].set(x[i]);
    }

    scalar error = run();
    if (this->log())
    {
        const unsigned int w = IOstream::defaultPrecision() + 7;
        this->logStream() << setf(ios_base::left);

        static bool writeHeader = true;
        if (writeHeader)
        {
            writeHeader = false;
            this->logStream()
                << setw(w) << "# x_0" << token::SPACE;
            for (label i = 1; i < x.size(); i++)
            {
                this->logStream()
                    << setw(w) << word("x_" + Foam::name(i))  << token::SPACE;
            }
            this->logStream() << setw(w) << "error" << endl;
        }
        forAll(x, i)
        {
            this->logStream()
                << setw(w) << x[i] << token::SPACE;
        }
        this->logStream() << setw(w) << error << endl;
    }

    return error;
}

