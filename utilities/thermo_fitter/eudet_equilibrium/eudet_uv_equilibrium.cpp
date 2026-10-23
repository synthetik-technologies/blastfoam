// eudet_uv_equilibrium: reference chemical equilibrium at fixed specific
// internal energy and density (UV problem) from the EuDetCode equilibrium
// solver, for the verification of the blastFoam equilibrium solver
// (tests/0D/verify_equilibrium.py).
//
// Usage:
//   eudet_uv_equilibrium <gasDb.json> <condensedDb.json> "<El:b ...>" \
//       rho1 e1 [rho2 e2 ...]
//
// The element abundances b are given in mol/kg (e.g. "N:54.75 O:14.57"), so
// that no ingredient lookup or default inert constraint (EuDetCode keeps 90%
// of an Al reactant inert) is involved; e is the absolute specific internal
// energy [J/kg] on the database's enthalpy datum. For each state the output is
//   state <rho> <e> <T> <p> <converged> <e of the solution>
//   <species> <mass fraction>      (all species)

#include <cstdlib>
#include <iomanip>
#include <map>
#include <iostream>
#include <sstream>
#include <string>

#include "ThermoChemistry/CombustionSolver.hpp"
#include "ThermoChemistry/CompositionUtils.hpp"

using namespace EUDetCode;

int main(int argc, char *argv[]) {
    if (argc < 6 || (argc - 4) % 2 != 0) {
        std::cerr << "Usage: " << argv[0]
                  << " gasDb condensedDb \"El:b ...\" rho e [rho e ...]\n";
        return 1;
    }

    // Element abundances
    std::vector<std::pair<std::string, double>> b;
    {
        std::stringstream ss(argv[3]);
        std::string tok;
        while (ss >> tok) {
            const auto colon = tok.find(':');
            b.emplace_back(tok.substr(0, colon),
                           std::atof(tok.substr(colon + 1).c_str()));
        }
    }

    // The reactant (one gas species per element) only sets up the solver
    // with the elements of the problem; the abundances are given directly
    const std::map<std::string, std::string> carrier = {
        {"O", "O2"}, {"N", "N2"}, {"C", "CO"}, {"H", "H2"}, {"Al", "AlO"}};
    std::string reactant;
    for (const auto &[el, v] : b) {
        reactant += carrier.at(el) + ":1 ";
    }
    CombustionSolver solver(reactant, "Ideal", defaultCompositionDbPaths(),
                            argv[1], argv[2], false);
    if (solver.hasAluminium()) {
        // No inert Al bound (EuDetCode default: 90% of the Al)
        solver.setInert("Al:0");
    }
    const MultiPhaseMixture &mix = solver.mixture();

    const std::vector<std::string> names = mix.elementNames();
    Eigen::VectorXd Aeq = Eigen::VectorXd::Zero(names.size());
    for (const auto &[el, v] : b) {
        const auto it = std::find(names.begin(), names.end(), el);
        if (it == names.end()) {
            std::cerr << "Element " << el << " not in the mixture\n";
            return 1;
        }
        Aeq[std::distance(names.begin(), it)] = v;
    }

    std::cout << std::setprecision(17);
    for (int i = 4; i < argc; i += 2) {
        const double rho = std::atof(argv[i]);
        const double e = std::atof(argv[i + 1]);
        CombustionSolver::CombustionResult res;
        try {
            res = solver.solveIsochoric(Aeq, 1.0 / rho, e);
        } catch (const std::exception &) {
            res.converged = false;
        }
        std::cout << "state " << rho << ' ' << e << ' ' << res.state.T << ' '
                  << res.state.P << ' ' << res.converged << ' '
                  << res.state.e << '\n';
        const Eigen::VectorXd &n = res.state.molalities;
        for (long k = 0; k < n.size(); ++k) {
            std::cout << mix.speciesName(k) << ' '
                      << n[k] * mix.speciesWeights()[k] << '\n';
        }
    }
    return 0;
}
