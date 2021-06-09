#include "dictionary.H"
#include "simpleEOS.H"
#include "OFstream.H"
#include "IFstream.H"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace Foam;

int main(int argc, char *argv[])
{
    fileName name("thermoDict");
    IFstream is(name);
    dictionary dict(is);

    //- Read state parameters

    scalar T = 0.0;

    autoPtr<simpleEOS> eosPtr(simpleEOS::New(dict));
    simpleEOS& eos = eosPtr();

    scalar p = 1e5;
    scalar rho = 1000;
    scalar e = 3.543e5;//eos.Es(rho, 0, T);
    scalar Cv = eos.Cv(rho, e, T);
    T = (e - eos.E(rho, e, T))/Cv;

    label n = 1000;
    std::vector<scalar> rhos(n, 0.0);
    std::vector<scalar> ps(n, 0.0);
    std::vector<scalar> Es(n, 0.0);
    forAll(rhos, i)
    {
        rhos[i] = 500 + 1500/scalar(n)*scalar(i);
        ps[i] = eos.p(rhos[i], e, T);
        Es[i] = eos.E(rhos[i], e, T);
    }

    e = eos.initializeEnergy(p, rho, e, T);
    Info<<"e: "<<eos.E(rho, e, T)<<" "<<e<<endl;
    Info<<"gamma: "<<eos.Gamma(rho, e, T)<<endl;
    Info<<"p: "<<eos.p(rho, e, T)<<endl;
    Info<<"c: "<<Foam::sqrt(eos.cSqr(p, rho, e, T))<<endl;
    Info<<"T: "<<eos.TRhoE(T,rho, e)<<endl;
    Info<<"rho: "<<eos.initializeRho(p, rho, e, T)<<endl;

//     plt::plot(rhos, Es, "g--");
//     plt::show();




    return 0;
}
