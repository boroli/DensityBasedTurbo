/* Reading the thermal properties for the solid region and calculating
   thermal coefficient. */

    Info<< "Reading thermal properties\n" << endl;

    IOdictionary thermalProperties
    (
        IOobject
        (
            "thermalProperties",
            runTime.constant(),
            stressMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Switch thermalStress(thermalProperties.lookup("thermalStress"));

    //  thermal coefficient
    dimensionedScalar threeKalpha
    (
        "threeKalpha",
        dimensionSet(0, 2, -2 , -1, 0),
        0
    );

    if (thermalStress)
    {
        // alpha=  coefficient of linear expansion
        dimensionedScalar alpha(thermalProperties.lookup("alpha"));

        Info<< "Calculating thermal coefficient\n" << endl;

        threeKalpha = threeK*alpha;

        Info<< "threeKalpha = " << threeKalpha.value() << " Pa/rhoS\n";
    }
