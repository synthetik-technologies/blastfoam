#include "dictionary.H"
#include "simpleBlastThermo.H"
#include "OFstream.H"
#include "IFstream.H"
#include "argList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList args(argc, argv, false, true);

    fileName name("thermoDict");
    IFstream is(name);
    dictionary dict(is);

    //- Read state parameters

    // scalar T = 60;

    autoPtr<simpleBlastThermo> eosPtr(simpleBlastThermo::New(dict));
    simpleBlastThermo& eos = eosPtr();

    scalar p = 1e5;
    scalar rho = 50;
    scalar T = 7.3;//eos.rhoPT(0.1, p, T);//1.225;
    // scalar e = eos.Es(rho, 0, T);

//     label n = 1000;
//     std::vector<scalar> rhos(n, 0.0);
//     std::vector<scalar> ps(n, 0.0);
//     std::vector<scalar> Es(n, 0.0);
//     forAll(rhos, i)
//     {
//         rhos[i] = 500 + 1500/scalar(n)*scalar(i);
//         ps[i] = eos.p(rhos[i], e, T);
//         Es[i] = eos.E(rhos[i], e, T);
//     }

    scalar e = eos.initializeEnergy(p, rho, e, T);
    // scalar T = eos.TRhoE(T, rho, e);

    Info<<"rho: "<< rho <<endl;
    Info<<"e: "<< e <<endl;
    Info<<"gamma: "<< eos.Gamma(rho, e, T) <<endl;
    Info<<"p: "<< eos.p(rho, e, T) <<endl;
    Info<<"c: "<< Foam::sqrt(eos.cSqr(p, rho, e, T)) <<endl;
    Info<<"T: "<< eos.TRhoE(T, rho, e) <<endl;
    Info<<"Cp: "<< eos.Cp(rho, e, T) <<endl;
    Info<<"Cv: "<< eos.Cv(rho, e, T) <<endl;
    Info<<"rho: "<< eos.rhoPT(rho, p, T) <<endl;
    Info<<"p: "<< eos.p(rho, e, T) <<endl;

    return 0;
}
