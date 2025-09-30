#include "fluidBlastThermo.H"
#include "cavitatingFluidBlastThermo.H"
#include "blendedBlastThermo.H"
#include "forDetBlastGases.H"
#include "makeDetBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

#include "constTransport.H"

#include "thermoModel.H"

#include "eConstBlastThermo.H"
#include "ePolynomialBlastThermo.H"

#include "idealGas.H"
#include "perfectGas.H"

#include "linearTillotson.H"

#include "specieBlast.H"
#include "rspecieBlast.H"

#include "forDetBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forCavvFluidEqns(uMu, rMu, uCp, rCp, uEos, Macro, Args...)            \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, idealGas, specieBlast, Macro, Args); \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, perfectGas, specieBlast, Macro, Args);


#define forCavlFluidEqns(uMu, rMu, uCp, rCp, Macro, Args...)                  \
    forCavvFluidEqns(uMu, rMu, uCp, rCp, linearTillotson, Macro, Args);

#define forCavvrFluidThermos(uMu, rMu, uCp, Macro, Args...)                    \
    forCavlFluidEqns(uMu, rMu, uCp, eConstThermo, Macro, Args);

#define forCavlFluidThermos(uMu, rMu, Macro, Args...)                         \
    forCavvrFluidThermos(uMu, rMu, eConstThermo, Macro, Args);

#define forCavvFluidTransports(uMu, Macro, Args...)                           \
    forCavlFluidThermos(uMu, constTransport, Macro, Args);

#define forCavlGasTransports(Macro, Args...)                                \
    forCavvFluidTransports(constTransport, Macro, Args);

#define forCavFluids(Macro, Args...)                                         \
    forCavlGasTransports(Macro, Args)

namespace Foam
{
    forCavFluids
    (
        makeDetThermo,
        fluidBlastThermo,
        cavitatingFluidBlastThermo,
        blendedBlastThermo
    );
}
