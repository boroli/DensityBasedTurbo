{
    Info << "Setting pressure" << endl;

    // interpolates the pressure from the fluid field to the
    // solid patches at the fluid-solid boundary
    scalarField solidPatchPressure =
        interpolatorFluidSolid.faceInterpolate
        (
            p.boundaryField()[fluidPatchID]
        );

    // tforce is Reference for traction displacement BC
    tForce.pressure() = solidPatchPressure;

    // calculates total force from pressure
    vector totalForcefromPressure =
        sum
        (
            p.boundaryField()[fluidPatchID]*
            mesh.Sf().boundaryField()[fluidPatchID]
        );

    Info << "Total force from pressure= " << totalForcefromPressure << endl;
}
