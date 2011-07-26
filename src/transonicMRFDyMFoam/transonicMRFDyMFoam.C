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

        runTime++;

        // rotate the mesh about the axis
        mesh.update();

#       include "solveFluid.H"

        runTime.write();

        Info<< "\n    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
