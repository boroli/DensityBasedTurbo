{
    volScalarField rhoSCs = rhoS*Cs;
    // solves temperature distribution  Tsolid with nonOrth correction
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
         solve
         (
                fvm::ddt(rhoSCs, Tsolid) - fvm::laplacian(Ks, Tsolid)
         );
    }
}
