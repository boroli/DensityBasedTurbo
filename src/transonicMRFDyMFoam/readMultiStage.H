    numberSubCycles = 1;
    numberRungeKuttaStages = 1;
    beta.setSize(numberRungeKuttaStages, 1.0);
    secondOrder = true;

if(mesh.solutionDict().found("multiStage"))
{
    // Read multi-stage coeffs dict
    dictionary multiStage = mesh.solutionDict().subDict("multiStage");

    // Read number of subCycles
    if (multiStage.found("numberSubCycles"))
    {
         multiStage.lookup("numberSubCycles") >> numberSubCycles;
    }

    // Read number of multi-stages
    if (multiStage.found("numberRungeKuttaStages"))
    {
         multiStage.lookup("numberRungeKuttaStages") >> numberRungeKuttaStages;
         beta.setSize(numberRungeKuttaStages, 1.0);
    }

    // Read multi-stage coeffs
    if(multiStage.found("RKCoeff") && numberRungeKuttaStages > 1)
    {
        ITstream rkCoeffInput = multiStage.lookup("RKCoeff");

        forAll(beta,betai)
        {
            rkCoeffInput >> beta[betai];
        }
    }
}

// read riemann solver coeffs
if(mesh.solutionDict().found("Riemann"))
{
    // Read Riemann coeffs dict
    dictionary riemann = mesh.solutionDict().subDict("Riemann");
    if (riemann.found("secondOrder"))
    {
        secondOrder = Switch(riemann.lookup("secondOrder"));
    }
}
