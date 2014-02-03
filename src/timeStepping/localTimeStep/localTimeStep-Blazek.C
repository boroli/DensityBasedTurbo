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

Class
    Foam::localTimeStep

Author
    Oliver Borm  All rights reserved.
    Sebastian Saegeler  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "localTimeStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::localTimeStep::localTimeStep
(
    const volVectorField& U,
    basicThermo& thermophysicalModel,
    compressible::turbulenceModel& turbulenceModel
)
:
    mesh_(U.mesh()),
    U_(U),
    thermophysicalModel_(thermophysicalModel),
    turbulenceModel_(turbulenceModel),
    CoDeltaT_
    (
        IOobject
        (
            "CoDeltaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("CoDeltaT", dimTime, 0.1),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaT_inv
    (
        IOobject
        (
            "deltaT_inv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaT_inv", dimTime,1.0),
        zeroGradientFvPatchVectorField::typeName
    ),
    deltaT_vis
    (
        IOobject
        (
            "deltaT_vis",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaT_vis", dimTime,1.0),
        zeroGradientFvPatchVectorField::typeName
   ),
    maxDeltaT_
    (
        IOobject
        (
            "maxDeltaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("maxDeltaT", dimTime,1.0),
        zeroGradientFvPatchVectorField::typeName
    ),
    maxDeltaTFluent_
    (
        IOobject
        (
            "maxDeltaTFluent",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("maxDeltaTFluent", dimTime,1.0),
        zeroGradientFvPatchVectorField::typeName
    )    
{};



void localTimeStep::update(scalar maxCo, Switch adjustTimeStep)
{

    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // interpolate velocity on cell faces
    surfaceVectorField Uf_
    (
        IOobject
        (
            "Uf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(U_))
    );
    // interpolate speed of sound on cell faces
    surfaceScalarField cf_
    (
        IOobject
        (
            "cf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(sqrt(thermophysicalModel_.Cp()/
        (thermophysicalModel_.Cv()*thermophysicalModel_.psi()))))
    );
    // interpolate density on cell faces
    surfaceScalarField rhof_
    (
        IOobject
        (
            "rhof",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(thermophysicalModel_.rho()))
    );
    // interpolate rho on cell faces
    surfaceScalarField muEfff_
    (
        IOobject
        (
            "muEfff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(turbulenceModel_.muEff()))
    );
    // interpolate kappa on cell faces
    surfaceScalarField kappaf_
    (
        IOobject
        (
            "kappaf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(thermophysicalModel_.Cp()/thermophysicalModel_.Cv()))
    );
    // minimum inviscid deltaT
    volScalarField deltaTMinValue_inv
    (
        IOobject
        (
            "deltaTMinValue_inv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaTMinValue_inv", dimTime,1.0),
        zeroGradientFvPatchVectorField::typeName
    );
    // minimum viscous deltaT
    volScalarField deltaTMinValue_vis
    (
        IOobject
        (
            "deltaTMinValue_vis",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaTMinValue_vis", dimTime,1.0),
        zeroGradientFvPatchVectorField::typeName
    );
    surfaceScalarField lambdaMax
    (
        IOobject
        (
            "lambdaMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(sqrt(thermophysicalModel_.Cp()/
        (thermophysicalModel_.Cv()*thermophysicalModel_.psi()))))
    );
    volScalarField phif_
    (
        IOobject
        (
            "phif",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("phif", dimensionSet(0,3,-1,0,0,0,0),SMALL),
        zeroGradientFvPatchVectorField::typeName
    );
    volScalarField LambdaVis_
    (
        IOobject
        (
            "LambdaVis",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("LambdaVis", dimensionSet(0,6,-1,0,0,0,0),SMALL),
        zeroGradientFvPatchVectorField::typeName
    );        
        
	volScalarField cellVolume
         (
            IOobject
            (
                "cellVolume",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimVolume,
            zeroGradientFvPatchScalarField::typeName
        );

    cellVolume.internalField() = mesh().V();
    cellVolume.correctBoundaryConditions();
	    
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();	    	   

    surfaceVectorField normal = Sf/magSf; //face normal vector

    // Inviscid time-step    
    forAll(Uf_, faceI)
    {
        label own = owner[faceI];
        //label nei = neighbour[faceI];

//	lambdaMax[faceI] = max(max((Uf_[faceI] & normal[faceI])+cf_[faceI],
//				(Uf_[faceI] & normal[faceI])-cf_[faceI]),(Uf_[faceI] & normal[faceI]));

//	phif_[own] += (magSf[faceI]*lambdaMax[faceI]);
	phif_[own] += (magSf[faceI]*(mag(Uf_[faceI]&normal[faceI])+cf_[faceI]));	     
    } 
    
    forAll(phif_.boundaryField(), patchI)
    {

        const fvPatchScalarField& pp = phif_.boundaryField()[patchI];
	
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];
	
	const fvsPatchVectorField& pUf = Uf_.boundaryField()[patchI];
	const fvsPatchScalarField& pcf = cf_.boundaryField()[patchI];
	
	fvsPatchScalarField& plambdaMax = lambdaMax.boundaryField()[patchI];	
	
	const unallocLabelList& faceCells =  pp.patch().faceCells();

        forAll(pp, facei)
        {
            vector normal = pSf[facei]/pMagSf[facei];

            label own = faceCells[facei];
    
//	    plambdaMax[facei]= max(max((pUf[facei] & normal)+pcf[facei],
//				(pUf[facei] & normal)-pcf[facei]),(pUf[facei] & normal));

//	    phif_[own] += (pMagSf[facei]*plambdaMax[facei]);					
	    phif_[own] += (pMagSf[facei]*(mag(pUf[facei]&normal)+pcf[facei]));
        }
    }       

//    deltaT_inv = cellVolume/phif_;


    // Viscous time-step
    if ( max(turbulenceModel_.muEff()).value() > SMALL )
    {    
	    forAll(rhof_, faceI)
	    {
	        label own = owner[faceI];
	        //label nei = neighbour[faceI];

		LambdaVis_[own] += sqr(magSf[faceI])*muEfff_[faceI]
					*(kappaf_[faceI]/rhof_[faceI]);	
     
    	    } 

	    forAll(LambdaVis_.boundaryField(), patchI)
	    {

	        const fvPatchScalarField& pp = LambdaVis_.boundaryField()[patchI];

	        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
	        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];
	
		const fvsPatchScalarField& pkappaf = kappaf_.boundaryField()[patchI];
		const fvsPatchScalarField& prhof = rhof_.boundaryField()[patchI];
		const fvsPatchScalarField& pmuEfff = muEfff_.boundaryField()[patchI];	
	
		const unallocLabelList& faceCells =  pp.patch().faceCells();

	        forAll(pp, facei)
	        {
	            vector normal = pSf[facei]/pMagSf[facei];

	            label own = faceCells[facei];
	    
	    	    LambdaVis_[own] += sqr(pMagSf[facei])
			*((pkappaf[facei]/prhof[facei])*pmuEfff[facei]);
        	}	
    	    }	
    }       

//    deltaT_vis = sqr(cellVolume)/LambdaVis_;
    
    volScalarField LambdaVIS = LambdaVis_/cellVolume; 

//    maxDeltaT_ = min(deltaT_vis, deltaT_inv);
    maxDeltaT_ = cellVolume/(phif_+4.0*LambdaVIS);
      
    CoDeltaT_ = maxCo* maxDeltaT_;


    if (adjustTimeStep)
    {
        CoDeltaT_ = min(CoDeltaT_);
    }

}


void localTimeStep::update(scalar maxCo)
{
    update(maxCo,"no");
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
