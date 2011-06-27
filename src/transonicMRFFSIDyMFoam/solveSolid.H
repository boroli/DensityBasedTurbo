
{
#   include "readStressedFoamControls.H"

    int iCorr = 0;
    scalar initialResidual = 0;

    do
    {
        volTensorField gradUsolid = fvc::grad(Usolid);

        if (thermalStress)
        {
#       include "solveTSolid.H"
        }

        // set-up equation of distortion Usolid [m]!!!
        fvVectorMatrix UEqn
        (
            fvm::d2dt2(Usolid)
         ==
            fvm::laplacian(2*mu + lambda, Usolid, "laplacian(DU,U)")
          + fvc::div
            (
                mu*gradUsolid.T() + lambda*(I*tr(gradUsolid)) - (mu + lambda)*gradUsolid,
                "div(sigma)"
            )
          - ( omega_ ^ (omega_ ^ stressMesh.C()))
        );

        // ads thermal expansion to UEqn
        if (thermalStress)
        {
            UEqn += threeKalpha*fvc::grad(Tsolid);
        }

        // solves equation of distortion Usolid
        initialResidual = UEqn.solve().initialResidual();

    } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

// calculates stress tensor
#       include "calculateStress.H"

}

