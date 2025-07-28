
volScalarField x
(
    IOobject
    (
        "x",
        runTime.timeName(),
        mesh
    ),
    mesh,
    0.0,
    "fixedValue"

);
volScalarField y
(
    IOobject
    (
        "y",
        runTime.timeName(),
        mesh
    ),
    mesh,
    0.0,
    "fixedValue"

);

volScalarField mixedx
(
    IOobject
    (
        "mixedx",
        runTime.timeName(),
        mesh
    ),
    mesh,
    0.0,
    "mixed"

);
volScalarField grady
(
    IOobject
    (
        "grady",
        runTime.timeName(),
        mesh
    ),
    mesh,
    0.0,
    "fixedGradient"

);

pointVectorField pointx
(
    IOobject
    (
        "pointx",
        runTime.timeName(),
        mesh
    ),
    pointMesh::New(mesh),
    0.5*vector::one,
    "fixedValue"

);
pointVectorField pointy
(
    IOobject
    (
        "pointy",
        runTime.timeName(),
        mesh
    ),
    pointMesh::New(mesh),
    3.0*vector::one,
    "fixedValue"
);
