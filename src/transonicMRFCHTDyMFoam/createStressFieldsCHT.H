    // Specific heat capacity for solid [J/kg/K]
    Info<< "Reading field Csolid\n" << endl;
    volScalarField Cs
    (
          IOobject
          (
              "C",
              runTime.timeName(),
              stressMesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE
          ),
          stressMesh
    );

    // Heat conductivity coefficient[kgm/s^3K] for solid mesh
    Info<< "Reading field Ksolid\n" << endl;
    volScalarField Ks
    (
          IOobject
          (
              "K",
              runTime.timeName(),
              stressMesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE
          ),
          stressMesh
    );

    //- Temperature [K]
    Info<< "Reading field Tsolid\n" << endl;
    volScalarField Tsolid
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            stressMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
         ),
         stressMesh
     );
