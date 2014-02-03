    label numberRungeKuttaStages = 1;
    scalarList beta(numberRungeKuttaStages, 1.0);
    bool secondOrder = true;
  

if(mesh.solutionDict().found("multiStage"))
{
    // Read multiStage Coeffs
    dictionary multiStage = mesh.solutionDict().subDict("multiStage");

    if (multiStage.found("numberRungeKuttaStages"))
    {
         multiStage.lookup("numberRungeKuttaStages") >> numberRungeKuttaStages;
         beta.setSize(numberRungeKuttaStages, 1.0);
    }

    if(multiStage.found("RKCoeff") && numberRungeKuttaStages > 1)
    {
        ITstream rkCoeffInput = multiStage.lookup("RKCoeff");

        forAll(beta,i)
        {
            rkCoeffInput >> beta[i];
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
