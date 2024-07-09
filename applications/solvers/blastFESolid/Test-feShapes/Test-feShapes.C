
#include "argList.H"
#include "Time.H"
#include "finiteElement.H"


using namespace Foam;


int main(int argc, char *argv[])
{
    argList::validArgs.append("transformations");
    argList::addOption("o", "Finite element order");
    argList::addOption("io", "Finite element order");

    #include "setRootCase.H"

    const word feType(args[1]);
    const finiteElement* elem = finiteElement::getRefFiniteElement
    (
        feType,
        args.optionLookupOrDefault<label>("o", 1)
    );
    const label io = args.optionLookupOrDefault<label>("io", 2);

    Info<<feType<<endl;
    const integrationRule& ir = elem->ir(io, io, io);
    forAll(ir, i)
    {
        Info<<"    ip: " << ir[i] << endl;
        Info<<"    shape: "<<elem->calcShape(elem->ir()[i])<<endl;
//         Info<<"    dshape: "<<elem->calcDShape(elem->ir()[i])<<endl;
    }

    return 0;
}
