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
    transonicSteadySRFFoam

Description
    Density-based compressible steady-state time-marching flow solver.

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"

#include "SRFModel.H"
#include "godunovFlux.H"
#include "localTimeStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
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

        Info<< "\n Time = " << runTime.value() << nl << endl;

        // Begin outer Loop for Runge - Kutta
        forAll (beta,i)
        {
            // Update primitive boundaries
            p.correctBoundaryConditions();
            h.correctBoundaryConditions();
            U.correctBoundaryConditions();

            // solve the approximate Riemann problem for this time step
            // reconstruct numerical fluxes at faces in a proper way
            Godunov.update(secondOrder);

            // get access to multi-stage coefficient in ddt scheme
            physDeltaT[0] = beta[i];

            // update local time step sizes just once each iteration for all ddt schemes
            localTimeStep.update(maxCo,adjustTimeStep);

            // \f$ \mu \left( \nabla \vec{U} + \left( \nabla \vec{U}
            // \right)^T - \frac{2}{nD} \left( \nabla \bullet \vec{U}
            // \right) I \right) \f$
            // nD = 3 in three dimensional space
            // for two- and one-dimensional computations this
            // implementation is wrong
            // is equal to -turbulence->devRhoReff(), but we do not need to
            // evaluate gradU twice
            volSymmTensorField tau
            (
                "tau",
                -turbulence->devRhoReff()
                -((2.0/3.0)*I)*rho*turbulence->k()
            );

            volScalarField rhoFlux = -fvc::div(Godunov.rhoFlux());

            volVectorField rhoUFlux = -fvc::div(Godunov.rhoUFlux())
                + fvc::div(tau);

            volScalarField rhoEFlux = -fvc::div(Godunov.rhoEFlux())
                // Viscous heating with
                + fvc::div( tau & U )
                // Fourier law with static enthalpy
                // with alphaEff() - the effective turbulent thermal diffusivity.
                // TODO: There seem to be sometimes problems with inviscid flows
                // when using harmonic interpolation
                + fvc::laplacian(turbulence->alphaEff(), h)
                // Molecular Diffusion and Turbulent Transport closure
                // Wilcox (2006): Chapter 5.4.3
                // should be better used DkEff(F1) instead of muEff(), but
                // this function is not virtual, now it is assumed that
                // \sigma_k = 5/3 is hard coded
                + fvc::laplacian
                  (
                      (turbulence->mu()+0.6*turbulence->mut()),
                      turbulence->k()
                  )
                ;

            rhoUFlux.internalField() -= rho.internalField()*SRF->Su()();
            rhoEFlux.internalField() -= rho.internalField()*(SRF->Fcentrifugal()() & U.internalField());

            // time integration
            solve
            (
                fvm::ddt(rho) == rhoFlux
            );

            // time integration
            solve
            (
                fvm::ddt(rhoU) == rhoUFlux
            );

            // time integration
            solve
            (
                fvm::ddt(rhoE) == rhoEFlux
            );

            // bound density
            boundMinMax(rho,rhoMin,rhoMax);

            // bound rhoE
            boundMinMax(rhoE,rhoEMin,rhoEMax);

            // Compute internal field of U
            U.dimensionedInternalField() =
                rhoU.dimensionedInternalField()
                /rho.dimensionedInternalField();

            // Update boundary conditions of U
            U.correctBoundaryConditions();

            // Bound the velocity
            volScalarField magU = mag(U);

            if (max(magU) > UMax)
            {
                Info<< "bounding " << U.name()
                    << " max: " << max(magU).value()
                    << endl;

                volScalarField Ulimiter = pos(magU - UMax)*UMax/(magU + smallU)
                    + neg(magU - UMax);
                Ulimiter.max(scalar(0));
                Ulimiter.min(scalar(1));

                U *= Ulimiter;
                U.correctBoundaryConditions();
            }

            const volScalarField kappa(thermo.Cp()/thermo.Cv());

            //Update static enthalpy:
            // The turbulent kinetic energy k is part of the total energy E
            // Therefore it needs to be subtracted from E in order to get
            // the static enthalpy h
            h = kappa*(rhoE/rho - 0.5*magSqr(U) - turbulence->k());

            // correct boundary conditions of static enthalpy
            h.correctBoundaryConditions();

            // bound enthalpy
            boundMinMax(h,hMin,hMax);

            // compute complete field of p
            p = (1.0 - 1.0/kappa)*rho*h;

            // correct boundary conditions of p
            p.correctBoundaryConditions();

            // bound pressure
            boundMinMax(p,pMin,pMax);

            // correct thermo physical properties
            // therefore one needs correct p and h fields
            thermo.correct();

            // Update boundary field of rho
            rho.boundaryField() = thermo.rho()().boundaryField();

            // Update boundary field of rhoU
            rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

            // Update boundary field of rhoE
            rhoE.boundaryField() =
                    rho.boundaryField()*
                    (
                        0.5*magSqr(U.boundaryField())
                      + turbulence->k()().boundaryField()
                    )
                  + p.boundaryField()/(kappa.boundaryField()-1.0);

            // needed for turbulence and CoEuler ddt scheme
            // and maybe elsewhere;
            phi = Godunov.rhoFlux();

            // Convergence check
            // Averaged L2-Norm of fluxes
            scalar L2NormRho = max(Foam::sqrt(gSum(sqr(rhoFlux.internalField())))
                /mesh.nCells(),SMALL);
            scalar L2NormRhoU = max(Foam::sqrt(gSum(magSqr(rhoUFlux.internalField())))
                    /mesh.nCells(),SMALL);
            scalar L2NormRhoE = max(Foam::sqrt(gSum(sqr(rhoEFlux.internalField())))
                /mesh.nCells(),SMALL);

            // Averaged L2-Norm of fluxes
            Info<< "rho residual: "
                << Foam::log10(L2NormRho)  << endl
                << "rhoU residual: "
                << Foam::log10(L2NormRhoU) << endl
                << "rhoE residual: "
                << Foam::log10(L2NormRhoE) << endl << endl;

        // End outer Loop for Runge - Kutta
        }

        //Update turbulence after the multi-stage time integration
        turbulence->correct();

        runTime.write();


        Info<< "\n    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
