
#   include "readMechanicalPropertiesCHT.H"

    // Youngs modulus [kg/(ms^2)]
    dimensionedScalar Em(mechanicalProperties.lookup("E"));
    // Poisson ratio [-]
    dimensionedScalar nuS(mechanicalProperties.lookup("nu"));

    Info<< "Normalising E: Em/rho\n" << endl;
    dimensionedScalar E = Em/rhoS;

    Info<< "Calculating Lame's coefficients\n" << endl;
    dimensionedScalar mu = E/(2.0*(1.0 + nuS));

    // calculating lambda and threeK for plane strain and 3D
    dimensionedScalar lambda = nuS*E/((1.0 + nuS)*(1.0 - 2.0*nuS));

    dimensionedScalar threeK = E/(1.0 - 2.0*nuS);

    // optional switch between planeStress and planeStrain
    Switch planeStress
    (
        mechanicalProperties.lookupOrDefault("planeStress",true)
    );

    // if planeStrees or planeStrain is selected in mechanicalProperties
    if (planeStress)
    {
        Info<< "Plane Stress\n" << endl;

        // Change lambda and threeK for plane stress
        lambda = nuS*E/((1.0 + nuS)*(1.0 - nuS));
        threeK = E/(1.0 - nuS);
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }

    Info << "mu = " << mu.value() << " " << mu.dimensions() << endl
         << "lambda = " << lambda.value() << " " << lambda.dimensions() << endl
         << "threeK = " << threeK.value() << " " << threeK.dimensions() << endl;

    // optionally read rotational velocity for additional centrifugal force
    dimensionedVector omega_
    (
        mechanicalProperties.lookupOrDefault
        (
            "omega",dimensionedVector("zeroVector",dimless/dimTime,vector::zero)
        )
    );

    Info << omega_ << endl;
