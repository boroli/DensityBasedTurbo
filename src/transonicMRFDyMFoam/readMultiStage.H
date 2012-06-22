
if(mesh.solutionDict().found("multiStage"))
{
    // Read multi-stage coeffs dict
    dictionary multiStage = mesh.solutionDict().subDict("multiStage");

    // Read number of subCycles
    if (multiStage.found("numberSubCycles"))
    {
         multiStage.lookup("numberSubCycles") >> numberSubCycles;
    }

    // Read multi-stage coeffs
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
