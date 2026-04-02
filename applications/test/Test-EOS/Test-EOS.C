#include "dictionary.H"
#include "simpleBlastThermo.H"
#include "OFstream.H"
#include "IFstream.H"
#include "argList.H"
// #include "inputSyntaxEntry.H"

using namespace Foam;

//- Convert keyword syntax to "dot" if the dictionary is "dot" syntax
word dotToSlash(const fileName& entryName)
{
    if
    (
       entryName.find('/') != string::npos
    )
    {
        wordList entryNames(entryName.components('/'));

        word entry(entryNames[0]);
        for (label i = 1; i < entryNames.size(); i++)
        {
            entry += word('.') + entryNames[i];
        }
        return entry;
    }
    else
    {
        return entryName;
    }
}


int main(int argc, char *argv[])
{
    argList::addOption("dict", "Dictionary to read from");
    argList::addOption("entry", "Path in dictionary to read from");
    argList::addBoolOption("calcDensity", "Calculate density");
    argList::addOption("specie", "Name of specie");
    argList::addOption("p", "Pressure [Pa]");
    argList::addOption("rho", "Density [kg/m^3]");
    argList::addOption("T", "Temperature [K]");
    argList args(argc, argv, false, true);

    const fileName dictPath
    (
        args.optionLookupOrDefault<fileName>("dict", "thermoDict")
    );
    fileName name(dictPath);
    IFstream is(name);
    dictionary parentDict(is);
    dictionary dict;

    word entryName;
    if (args.optionReadIfPresent("entry", entryName))
    {
        const entry* entryPtr = parentDict.lookupScopedEntryPtr
        (
            entryName,
            false,
            true            // Support wildcards
        );
        if (!entryPtr)
        {
            wordList cmpts(fileName(entryName).components('/'));
            word scopedName;
            while (cmpts.size() && !entryPtr)
            {
                cmpts.setSize(cmpts.size()-1);
                scopedName = cmpts[0];
                for (label i = 1; i < cmpts.size(); i++)
                {
                    scopedName += word('.') + cmpts[i];
                }
                entryPtr = parentDict.lookupScopedEntryPtr
                (
                    scopedName,
                    false,
                    true            // Support wildcards
                );
            }
            FatalErrorInFunction
                << "Could not find " << scopedName << endl
                << abort(FatalError);
        }
        dict = entryPtr->dict();
    }
    else
    {
        dict = parentDict;
    }
    //- Read state parameters

    // scalar T = 60;

    autoPtr<simpleBlastThermo> eosPtr
    (
        simpleBlastThermo::New
        (
            dict,
            args.optionLookupOrDefault("specie", word::null)
        )
    );
    simpleBlastThermo& eos = eosPtr();

    scalar p = -1;
    scalar rho = -1;
    scalar T = -1;
    scalar e = -1;
    if (args.optionFound("p")) p = args.optionRead<scalar>("p");
    else dict.readIfPresent("p", p);

    if (args.optionFound("rho")) rho = args.optionRead<scalar>("rho");
    else dict.readIfPresent("rho", rho);

    if (args.optionFound("T")) T = args.optionRead<scalar>("T");
    else dict.readIfPresent("T", T);

    const label nRead = label(p > 0) + label(rho > 0) + label(T > 0);
    if (nRead == 3)
    {
        dict.lookupOrDefault("calculateDensity", false);
        rho = eos.rhoPT(rho, p, T);
    }
    else if (nRead < 2)
    {
        FatalErrorInFunction
            << nRead<<" Atleast 2 of p, rho, T must be specified" << endl
            << abort(FatalError);
    }
    else if (rho < 0)
    {
        rho = eos.rhoPT(rho, p, T);
        e = eos.initializeEnergy(p, rho, e, T);
    }
    else if (p < 0)
    {
        e = eos.Es(rho, 0, T);
        p = eos.p(rho, e, T);
    }
    else if (T < 0)
    {
        T = 300.0;
        e = eos.initializeEnergy(p, rho, e, T);
        T = eos.TRhoE(T, rho, e);
    }

    Info<< "Initial values: " << nl << incrIndent
        << indent << "p: " << p << nl
        << indent << "e: " << e << nl
        << indent << "rho: " << rho << nl
        << indent << "T: " << T << nl
        << endl << decrIndent;

    Info<< "Derived and recalculated: " << nl << incrIndent
        << indent <<"gamma: "<< eos.Gamma(rho, e, T) + 1.0 << nl
        << indent <<"speed of sound: "<< Foam::sqrt(eos.cSqr(p, rho, e, T)) << nl
        << indent <<"Cp: "<< eos.Cp(rho, e, T) << nl
        << indent <<"Cv: "<< eos.Cv(rho, e, T) << nl
        << indent <<"dpdT: "<< eos.dpdT(rho, e, T) << nl
        << indent <<"dpdv: "<< eos.dpdv(rho, e, T) << nl
        << indent << "p(rho, T): " << eos.p(rho, e, T) << nl
        << indent << "Es(rho, T): " << eos.Es(rho, e, T) << nl
        << indent << "e(p, rho, T): " << eos.initializeEnergy(p, rho, e, T) << nl
        << indent << "rho(p, T): " << eos.rhoPT(rho, p, T) << nl
        << indent << "T(rho, e): " << eos.TRhoE(T, rho, e) << nl
        << indent << "mu(rho, e, T): " << eos.mu(rho, e, T) << nl
        << indent << "kappa(rho, e, T): " << eos.kappa(rho, e, T) << nl
        << indent << "Pr(rho, e, T): " << eos.mu(rho, e, T)*eos.Cp(rho, e, T)/eos.kappa(rho, e, T) << nl
        << endl;

    return 0;
}
