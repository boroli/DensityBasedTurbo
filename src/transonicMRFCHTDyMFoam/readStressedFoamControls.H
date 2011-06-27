    const dictionary& stressControl =
        stressMesh.solutionDict().subDict("stressedFoam");

    // number of corrector loops of calculation of Usolid
    int nCorr(stressControl.lookupOrDefault("nCorrectors",2));

    // number of corrector loops for non orthogonal errors
    // of temperature distribution
    int nNonOrthCorr(stressControl.lookupOrDefault("nNonOrthogonalCorrectors",2));

    // convergence tolerance of calculation of Usolid
    scalar convergenceTolerance(readScalar(stressControl.lookup("U")));
