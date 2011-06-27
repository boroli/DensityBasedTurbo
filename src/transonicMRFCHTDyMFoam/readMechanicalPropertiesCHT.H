   Info<< "Reading mechanical properties\n" << endl;

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            stressMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Density [kg/m^3]
    dimensionedScalar rhoS(mechanicalProperties.lookup("rho"));

