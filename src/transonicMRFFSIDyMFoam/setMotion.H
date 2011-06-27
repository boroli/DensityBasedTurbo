{
// Setting mesh motion

    // takes the data on the solid patch (cell centres)
    // and interpolates it into points
    pointVectorField solidPointsDispl =
        cpi.interpolate(Usolid - Usolid.oldTime());

    word ptDisplacementName("volPointInterpolate((U-U_0))");

    vectorField newPoints =
        stressMesh.points()
      + solidPointsDispl.internalField();

    Info << "Move stress points" << endl;
    stressMesh.movePoints(newPoints);

//     Info << "Move stress points" << endl;
//     // use solidDisplacement motion solver for point motion
//     word ptDisplacementName("pointDisplacement");
//     stressMesh.update();

//++++++++++++++++++++++++++++++++++//
//     interpolatorSolidFluid.movePoints();
//++++++++++++++++++++++++++++++++++//

    const pointVectorField& ptDisplacement =
        stressMesh.lookupObject<pointVectorField>(ptDisplacementName);

    //takes the point data and interpolates from solid to fluid
    vectorField fluidPatchPointsDispl =
        interpolatorSolidFluid.pointInterpolate
        (
            ptDisplacement.boundaryField()[solidPatchID].
            patchInternalField()
        );
//++++++++++++++++++++++++++++++++++//
// laplaceFaceDecomposition motion solver

//     // takes the point data (that is mesh points) and interpolates into 
//     // whatever the motion solver wants: for cell decomposition it just uses cell data,
//     // for face decomposition it will do a face centre as well in a consistent manner
//     motionUFluidPatch ==
//         tppi.pointToPointInterpolate
//         (
//             fluidPatchPointsDispl/runTime.deltaT().value()
//         );

//++++++++++++++++++++++++++++++++++//
// Mesquite motion solver

    // Boundary motion specified for the motionSolver
    pointField& refPoints = const_cast<pointField&>
    (
        mesh.lookupObject<pointField>("refPoints")
    );

    // Assign boundary conditions to the motion solver at patch 'pI'
    const labelList& meshPts = mesh.boundaryMesh()[fluidPatchID].meshPoints();

    // Apply displacement
    forAll(meshPts,pointI)
    {
        refPoints[meshPts[pointI]] += fluidPatchPointsDispl[pointI];
    }
//++++++++++++++++++++++++++++++++++//

    Info << "Move fluid points" << endl;
    mesh.update();

// calculates volume continuity errors
#   include "volContinuity.H"

    Info << "Motion magnitude: mean = "
        << average(mag(Usolid.boundaryField()[solidPatchID]))
        << " max = " 
        << max(mag(Usolid.boundaryField()[solidPatchID])) << endl;
}
