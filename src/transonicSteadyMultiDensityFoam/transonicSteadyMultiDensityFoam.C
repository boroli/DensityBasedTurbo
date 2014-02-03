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
    transonicSteadyMultiDensityFoam

Description
    Density-based compressible time-marching flow solver for Roe-Riemannsolver.
    Preconditioning-matrix and implementation according to
    J. Weiss & Wayne Smith, AIAA Journal, Vol.33, No.11, 1995; and
    J. Blazek, Computational Fluid Dynamics: Principles and Applications, 2nd Ed, 2006
    

Author
    Sebastian Saegeler  All rights reserved.
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
//#include "RASModel.H"
#include "turbulenceModel.H"

#include "scalarMatrices.H"

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
#       include "readMultiStage.H"
#       include "readFieldBounds.H"

        runTime++;
	
        // Update local time step sizes just once each iteration for all ddt schemes.
	// Local time-stepping-scheme should be consistent with determination of the
	// characteristic lenght (deltaX), which is necessary for viscous velocity
        localTimeStep.update(maxCo,adjustTimeStep);
	volScalarField deltaX = localTimeStep.deltaX();	

        Info<< "\n Time = " << runTime.value() << nl << endl;

	// Begin outer Loop for Runge - Kutta
    	forAll (beta,i)
    	{	
	
    		// Update primitive boundaries
    		p.correctBoundaryConditions();
		U.correctBoundaryConditions();
    		T.correctBoundaryConditions();

		// solve the approximate Riemann problem for this time step
		// reconstruct numerical fluxes at faces in a proper way
        	Godunov.update(secondOrder);
			
                // get access to multi-stage coefficient in ddt scheme
                physDeltaT[0] = beta[i];
		
                // reynolds stress
            	volSymmTensorField tau
            	(
                    "tau",
                    -turbulence->devRhoReff()
                    -((2.0/3.0)*I)*rho*turbulence->k()
            	);				
	

		volScalarField rhoFlux  = fvc::div(Godunov.rhoFlux());
		
		volVectorField rhoUFlux = fvc::div(Godunov.rhoUFlux())
						- fvc::div(tau);
						
		volScalarField rhoEFlux = fvc::div(Godunov.rhoEFlux())
                				// Viscous heating with
                				- fvc::div( tau & U )
                				// Fourier law with static enthalpy
                				// with alphaEff() - the effective turbulent thermal diffusivity.
                				// TODO: There seem to be sometimes problems with inviscid flows
                				// when using harmonic interpolation
                				- fvc::laplacian(turbulence->alphaEff(), h)
                				// Molecular Diffusion and Turbulent Transport closure
                				// Wilcox (2006): Chapter 5.4.3
                				// should be better used DkEff(F1) instead of muEff(), but
                				// this function is not virtual, now it is assumed that
                				// \sigma_k = 1.0 is hard coded
                				- fvc::laplacian
                  				(
                      					(turbulence->mu()+1.0*turbulence->mut()),  // turbulence->DkEff()
                      					turbulence->k()
                  				)				
                				;

		// Evaluate GammaInverse at each cell
		// perform multiplication of inverse Jacobian with conservative fluxes
		volScalarField Cp = thermo.Cp();
		volScalarField Cv = thermo.Cv();
		volScalarField R  = Cp-Cv;
		kappa = max(Cp/Cv,SMALL);

		// speed of sound
		c = sqrt(kappa*R*T);

		// Calculate the total Enthalpy		
		volScalarField H = h + 0.5*magSqr(U);//+ turbulence->k();
		H.correctBoundaryConditions();

		// deviation of density by temperature
		volScalarField DrhoDT = -rho/T;
		
		// - find the viscous velocity
		volScalarField deltaX = localTimeStep.deltaX();
		volScalarField visVel = turbulence->muEff()/(deltaX*rho);
			
		// loop over all cells
		forAll(rhoFlux, celli)
		{	
			// create a 5x5 matrix
			scalarRectangularMatrix GammaInverse(5,5);
			
			// reference velocity
//			Ur[celli] = min(max(mag(U[celli]), 0.0001*c[celli]), c[celli]);
			Ur[celli] = min(max(max(mag(U[celli]), 0.1*c[celli]),visVel[celli]), c[celli]);
//			Ur[celli] = min(max(max(mag(U[celli]), 0.0001*c[celli]),visVel[celli]), c[celli]);			
//			Ur[celli] = min(max(max(max(visVel[celli],mag(U[celli])), SMALL*c[celli]), heatVel[celli]), c[celli]);

			scalar Teta = scalar(1.0)/sqr(Ur[celli]) + scalar(1.0)/(T[celli]*Cp[celli]);
			
			const scalar rRho = scalar(1.0)/rho[celli];
			
			scalar Psi = 1.0/(Teta*sqr(c[celli])-(kappa[celli]-1.0));
			scalar a3  = (kappa[celli]-1.0)*Psi*T[celli]*rRho;
			scalar a2  = (kappa[celli]-1.0)*Psi;

			// filling up the inverse of Gamma
			GammaInverse[0][0] =  a2* (sqr(c[celli])/(kappa[celli]-1)-(H[celli]-magSqr(U[celli])));
			GammaInverse[0][1] = -a2* U[celli][0];
			GammaInverse[0][2] = -a2* U[celli][1];
			GammaInverse[0][3] = -a2* U[celli][2];
			GammaInverse[0][4] =  a2;
		
			GammaInverse[1][0] = -U[celli][0]*rRho;
			GammaInverse[1][1] =  rRho;
		
			GammaInverse[2][0] = -U[celli][1]*rRho;
			GammaInverse[2][2] =  rRho;
		
			GammaInverse[3][0] = -U[celli][2]*rRho;
			GammaInverse[3][3] =  rRho;
		
			GammaInverse[4][0] =  a3* (1.0-Teta*(H[celli]-magSqr(U[celli])));
			GammaInverse[4][1] = -Teta* a3* U[celli][0];
			GammaInverse[4][2] = -Teta* a3* U[celli][1];
			GammaInverse[4][3] = -Teta* a3* U[celli][2];
			GammaInverse[4][4] =  Teta* a3;
			
		
			// perform matrix vector multiplication
 			pResidual[celli]    = GammaInverse[0][0]*rhoFlux[celli] + GammaInverse[0][1]*rhoUFlux[celli][0] + GammaInverse[0][2]*rhoUFlux[celli][1] + GammaInverse[0][3]*rhoUFlux[celli][2] + GammaInverse[0][4]*rhoEFlux[celli];
 	
 			UResidual[celli][0] = GammaInverse[1][0]*rhoFlux[celli] + GammaInverse[1][1]*rhoUFlux[celli][0];
 		
 			UResidual[celli][1] = GammaInverse[2][0]*rhoFlux[celli] + GammaInverse[2][2]*rhoUFlux[celli][1];
 		
 			UResidual[celli][2] = GammaInverse[3][0]*rhoFlux[celli] + GammaInverse[3][3]*rhoUFlux[celli][2];
 	
 			TResidual[celli]    = GammaInverse[4][0]*rhoFlux[celli] + GammaInverse[4][1]*rhoUFlux[celli][0] + GammaInverse[4][2]*rhoUFlux[celli][1] + GammaInverse[4][3]*rhoUFlux[celli][2] + GammaInverse[4][4]*rhoEFlux[celli];
		
		}
	
 		pResidual.correctBoundaryConditions();
 		UResidual.correctBoundaryConditions();
 		TResidual.correctBoundaryConditions();
	
		// Update boundary field and values
    		forAll(p.boundaryField(), patchi)
    		{
			fvPatchScalarField& pPResidual = pResidual.boundaryField()[patchi];
			fvPatchVectorField& pUResidual = UResidual.boundaryField()[patchi];
			fvPatchScalarField& pTResidual = TResidual.boundaryField()[patchi];

			const fvPatchScalarField& prho = rho.boundaryField()[patchi];
	
// 			const fvPatchScalarField& pp = p.boundaryField()[patchi];
			const fvPatchVectorField& pU = U.boundaryField()[patchi];
 			const fvPatchScalarField& pT = T.boundaryField()[patchi];
	
			const fvPatchScalarField& pRhoFlux  = rhoFlux.boundaryField()[patchi];
			const fvPatchVectorField& pRhoUFlux = rhoUFlux.boundaryField()[patchi];
			const fvPatchScalarField& pRhoEFlux = rhoEFlux.boundaryField()[patchi];
	
			const fvPatchScalarField& pCp = Cp.boundaryField()[patchi];
			const fvPatchScalarField& pH  = H.boundaryField()[patchi];
			const fvPatchScalarField& pc  = c.boundaryField()[patchi];
			const fvPatchScalarField& pkappa = kappa.boundaryField()[patchi];
	
//			const fvPatchScalarField& pDrhoDT = DrhoDT.boundaryField()[patchi];			
//			const fvPatchScalarField& pTeta   = Teta.boundaryField()[patchi];
			fvPatchScalarField& pUr	  = Ur.boundaryField()[patchi];
			const fvPatchScalarField& pvisVel  = visVel.boundaryField()[patchi];
//			const fvPatchScalarField& pheatVel = heatVel.boundaryField()[patchi];
			

 			forAll(prho, facei)
         		{
							
				scalarRectangularMatrix GammaInverse(5,5);
				
//				pUr[facei] = min(max(mag(pU[facei]), 0.0001*pc[facei]), pc[facei]);
				pUr[facei] = min(max(max(mag(pU[facei]), 0.1*pc[facei]),pvisVel[facei]), pc[facei]);
//				pUr[facei] = min(max(max(mag(pU[facei]), 0.0001*pc[facei]),pvisVel[facei]), pc[facei]);
//				pUr[facei] = min(max(max(max(pvisVel[facei],mag(pU[facei])), SMALL*pc[facei]), pheatVel[facei]), pc[facei]);				
			
				scalar Teta = scalar(1.0)/sqr(pUr[facei]) + scalar(1.0)/(pT[facei]*pCp[facei]);
				
				scalar const rRho = scalar(1.0)/prho[facei];				
			
				scalar Psi = 1.0/(Teta*sqr(pc[facei])-(pkappa[facei]-1.0));
				scalar a3  = (pkappa[facei]-1.0)*Psi*pT[facei]*rRho;
				scalar a2  = (pkappa[facei]-1.0)*Psi;

				// filling up the inverse of Gamma at the boundaries
				GammaInverse[0][0] =  a2* (sqr(pc[facei])/(pkappa[facei]-1)-(pH[facei]-magSqr(pU[facei])));
				GammaInverse[0][1] = -a2* pU[facei][0];
				GammaInverse[0][2] = -a2* pU[facei][1];
				GammaInverse[0][3] = -a2* pU[facei][2];
				GammaInverse[0][4] =  a2;
		
				GammaInverse[1][0] = -pU[facei][0]*rRho;
				GammaInverse[1][1] =  rRho;
		
				GammaInverse[2][0] = -pU[facei][1]*rRho;
				GammaInverse[2][2] =  rRho;
		
				GammaInverse[3][0] = -pU[facei][2]*rRho;
				GammaInverse[3][3] =  rRho;
		
				GammaInverse[4][0] =  a3* (1.0-Teta*(pH[facei]-magSqr(pU[facei])));
				GammaInverse[4][1] = -Teta* a3* pU[facei][0];
				GammaInverse[4][2] = -Teta* a3* pU[facei][1];
				GammaInverse[4][3] = -Teta* a3* pU[facei][2];
				GammaInverse[4][4] =  Teta* a3;				
				

				// perform matrix vector multiplication
	 			pPResidual[facei]    = GammaInverse[0][0]*pRhoFlux[facei] + GammaInverse[0][1]*pRhoUFlux[facei][0] + GammaInverse[0][2]*pRhoUFlux[facei][1] + GammaInverse[0][3]*pRhoUFlux[facei][2] + GammaInverse[0][4]*pRhoEFlux[facei];
 	
 				pUResidual[facei][0] = GammaInverse[1][0]*pRhoFlux[facei] + GammaInverse[1][1]*pRhoUFlux[facei][0];
 			
 				pUResidual[facei][1] = GammaInverse[2][0]*pRhoFlux[facei] + GammaInverse[2][2]*pRhoUFlux[facei][1];
 			
 				pUResidual[facei][2] = GammaInverse[3][0]*pRhoFlux[facei] + GammaInverse[3][3]*pRhoUFlux[facei][2];
 		
 				pTResidual[facei]    = GammaInverse[4][0]*pRhoFlux[facei] + GammaInverse[4][1]*pRhoUFlux[facei][0] + GammaInverse[4][2]*pRhoUFlux[facei][1] + GammaInverse[4][3]*pRhoUFlux[facei][2] + GammaInverse[4][4]*pRhoEFlux[facei];
 			}			
    		}		
	
		// time integration of p		
        	solve
        	(
            	    fvm::ddt(p) == -pResidual
        	);

		// time integration of U
        	solve
        	(
            	    fvm::ddt(U) == -UResidual
        	);
	
		// time integration of T
        	solve
        	(
            	    fvm::ddt(T) == -TResidual
        	);
//        	Info<< "End Solving Equationset" << endl;

		// not sure if this is needed here		
    		p.correctBoundaryConditions();
		U.correctBoundaryConditions();
    		T.correctBoundaryConditions();

/*****************************************************************/
//	This can be commented out if using a Tpsi thermo class,
//	otherwhise T-field is evaluated twice

		//Update static enthalpy, only valid for perfect gases
		h = thermo.Cp()*T;
    		// correct boundary conditions of static enthalpy
    		h.correctBoundaryConditions();
		
    		// bound enthalpy
    		//boundMinMax(h,hMin,hMax);		
		
/****************************************************************/
		
    		// correct thermo physical properties
    		thermo.correct();
			
    		// update conservative variables
    		rho = thermo.rho();
		
    		rhoU = rho*U;
    		rhoE = rho*(h + 0.5*magSqr(U))-p; //+ rho*turbulence->k());
    
    		// needed for turbulence, and CoEuler time integration
    		phi = Godunov.rhoFlux();
		
		// Convergence check
                // Averaged L2-Norm of fluxes
		scalar L2NormRho = max(Foam::sqrt(sum(sqr(rhoFlux.internalField()))
                			/mesh.nCells()),SMALL);		
                scalar L2NormP = max(Foam::sqrt(sum(sqr(pResidual.internalField()))
                			/mesh.nCells()),SMALL);
                scalar L2NormU = max(Foam::sqrt(sum(magSqr(UResidual.internalField()))
                			/mesh.nCells()),SMALL);
                scalar L2NormT = max(Foam::sqrt(sum(sqr(TResidual.internalField()))
                			/mesh.nCells()),SMALL);

            // Averaged L2-Norm of fluxes
            Info<< "rho residual: "
                << Foam::log10(L2NormRho)  << endl
	    	<< "p residual: "
                << Foam::log10(L2NormP)  << endl
                << "U residual: "
                << Foam::log10(L2NormU) << endl
                << "T residual: "
                << Foam::log10(L2NormT) << endl << endl;
			
    	} // End outer Loop for Runge - Kutta

        //Update turbulence after the multi-stage time integration
        turbulence->correct();
	
	//Some monitor-variables
	dimensionedScalar VMax  = max(mag(U));
	dimensionedScalar MaMax = max(mag(U)/c);
	dimensionedScalar TtMax = max(T + 0.5*magSqr(U)/thermo.Cp());// + turbulence->k())/thermo.Cp());
	dimensionedScalar ptMax = max(p + 0.5*rho*magSqr(U));
	
	Info<< "\n VMax = " << VMax.value() << "  MaMax = " << MaMax.value() << "  TtMax = " << TtMax.value() << "  ptMax = " << ptMax.value() << endl << endl;

        
	// Output
        runTime.write();	

        Info<< "\n    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
