/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    transonicMRFDyMFoam

Description
    Density-based compressible multi-stage (cf. JST) time-marching flow solver.
    The multi-stage coefficients and the number of multi-stages are runtime
    selectable. In order to evaluate the convective terms, a Riemann solver is
    utilized.

    MRF and dynamic prescribed mesh motion is also implemented.

    TODO: Make the Riemann solver and multidimensional limiter run-time selectable.

    References:
    Jameson, Schmidt and Turkel, "Numerical Solution of the Euler Equations by
    Finite Volume Methods Using Runge-Kutta Time-Stepping Schemes",
    AIAA 1981-1259

    This is a special version of the unsteady solver and does only work, if the
    mesh motion is a rigid body rotation about the global z-axis! You have been
    warned. In that case the velocity components are stored before the discrete
    mesh motion as radial and circumferential components. After the rotation
    the cartesian components are computed at this new position. This is done,
    because now the time derivative in the momentum equation can now be
    considered as it is in the relative frame of reference, as the vector did
    not change during the mesh motion according to the relative frame of
    reference. Thus the additional source term ( \omega \times (\rho \vec{U}) )
    can be applied. This should speed up the computation, as one will need less
    inner iterations.

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"

#include "MRFZones.H"
#include "dynamicFvMesh.H"
#include "godunovFlux.H"

#include "localTimeStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"

     label rhoUOldFieldI = 0;
     volVectorField& rhoUOld    = vectorFieldList[rhoUOldFieldI];
     volVectorField& rhoUOldOld = vectorFieldList[rhoUOldFieldI+vectorFieldNames.size()];

    const volScalarField& rhoCrOld =
        mesh.lookupObject<volScalarField>("rhoCrOld");
    const volScalarField& rhoCrOldOld =
        mesh.lookupObject<volScalarField>("rhoCrOldOld");

    const volScalarField& rhoCuOld =
        mesh.lookupObject<volScalarField>("rhoCuOld");
    const volScalarField& rhoCuOldOld =
        mesh.lookupObject<volScalarField>("rhoCuOldOld");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"

// this definition of compressibleCourantNo is not correct
// #       include "compressibleCourantNo.H"

// adjusting of physical deltaT is no longer needed with local-time stepping
// #       include "setDeltaT.H"

#       include "readMultiStage.H"
#       include "readFieldBounds.H"

        if (adjustTimeStep)
        {
            localTimeStep.update(maxCo,adjustTimeStep);
            runTime.setDeltaT
            (
                min
                (
                    min(localTimeStep.CoDeltaT()).value(),
                    maxDeltaT
                )
            );
            numberSubCycles = 1;
        }

        runTime++;

        // rotate the mesh about the axis
        mesh.update();

#       include "solveUnsteadyFluid.H"

        runTime.write();

        Info<< "\n    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
