
# include "createStressFieldsCHT.H"

    // Displacement [m]
    Info<< "Reading field U\n" << endl;
    volVectorField Usolid
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            stressMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        stressMesh
    );

    // interpolates volume fields to point fields for moving mesh
    const volPointInterpolation& cpi = volPointInterpolation::New(stressMesh);
