// Initialise solid field pointer lists

// General fields
PtrList<solidModel> solidModels(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh properties for region "
    << solidRegions[i].name() << nl << endl;

    solidModels.set
    (
        i,
        solidModel::New(solidRegions[i])
    );
}
