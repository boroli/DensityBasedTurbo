    label numberRungeKuttaStages = 1;
    scalarList beta(1, 1.0);
    bool secondOrder = true;

if(mesh.solutionDict().found("multiStage"))
{
    // Read multiStage Coeffs
    dictionary multiStage = mesh.solutionDict().subDict("multiStage");

    if(multiStage.found("RKCoeff"))
    {
         multiStage.lookup("RKCoeff") >> beta;
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
