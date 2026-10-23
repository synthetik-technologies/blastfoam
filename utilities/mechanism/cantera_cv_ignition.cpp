// cantera_cv_ignition: constant-volume, adiabatic ignition of an ideal gas
// with a Cantera mechanism (reference for tests/0D/verify_gas_kinetics.py).
//
// Usage:
//   cantera_cv_ignition <mechanism.yaml> <T0> <p0> "<Sp:Y, ...>" <tEnd> [nOut]
//
// Output: "time T p" rows (nOut + 1 equidistant times and every internal
// step), then "ignition <time of max dT/dt>", "final T p" and the final mass
// fractions "Y <species> <value>".

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "cantera/core.h"
#include "cantera/zerodim.h"

using namespace Cantera;

int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0]
                  << " mech.yaml T0 p0 \"Sp:Y, ...\" tEnd [nOut]\n";
        return 1;
    }
    auto sol = newSolution(argv[1], "", "none");
    auto gas = sol->thermo();
    gas->setState_TPY(std::atof(argv[2]), std::atof(argv[3]), argv[4]);
    const double tEnd = std::atof(argv[5]);
    const int nOut = argc > 6 ? std::atoi(argv[6]) : 200;

    IdealGasReactor r;
    r.insert(sol);
    ReactorNet net;
    net.addReactor(r);
    net.setTolerances(1e-10, 1e-16);

    std::cout << std::setprecision(10);
    double tPrev = 0, TPrev = r.temperature(), dTdtMax = 0, tIgn = 0;
    double tNext = tEnd/nOut;
    std::cout << 0.0 << ' ' << r.temperature() << ' ' << r.pressure() << '\n';
    while (net.time() < tEnd) {
        const double t = net.step();
        const double T = r.temperature();
        const double dTdt = (T - TPrev)/(t - tPrev);
        if (dTdt > dTdtMax) {
            dTdtMax = dTdt;
            tIgn = 0.5*(t + tPrev);
        }
        tPrev = t;
        TPrev = T;
        if (t >= tNext || t >= tEnd) {
            std::cout << t << ' ' << T << ' ' << r.pressure() << '\n';
            tNext += tEnd/nOut;
        }
    }
    std::cout << "ignition " << tIgn << '\n';
    std::cout << "final " << r.temperature() << ' ' << r.pressure() << '\n';
    const auto& th = r.contents();
    for (size_t k = 0; k < th.nSpecies(); k++) {
        std::cout << "Y " << th.speciesName(k) << ' ' << th.massFraction(k) << '\n';
    }
    return 0;
}
