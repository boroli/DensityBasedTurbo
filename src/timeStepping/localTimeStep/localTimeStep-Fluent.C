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
//    deltaX_
//    (
//        IOobject
//        (
//            "deltaX",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("deltaX", dimLength, SMALL),
//        zeroGradientFvPatchScalarField::typeName
//    ),
//    Lambda_inv
//    (
//        IOobject
//        (
//            "Lambda_inv",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("Lambda_inv", dimensionSet(0,3,1,0,0,0,0),1.0),
//        zeroGradientFvPatchVectorField::typeName
//    ),
//    Lambda_vis
//    (
//        IOobject
//        (
//            "Lambda_vis",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("Lambda_vis", dimensionSet(0,3,1,0,0,0,0),1.0),
//        zeroGradientFvPatchVectorField::typeName
//    ),
//    deltaT_inv
//    (
//        IOobject
//        (
//            "deltaT_inv",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("deltaT_inv", dimTime,1.0),
//        zeroGradientFvPatchVectorField::typeName
//    ),
//    deltaT_vis
//    (
//        IOobject
//        (
//            "deltaT_vis",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("deltaT_vis", dimTime,1.0),
//        zeroGradientFvPatchVectorField::typeName
//    ),
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
//{
//    updateDeltaX();
//};

//void localTimeStep::updateDeltaX()
//{
//    const unallocLabelList& owner = mesh().owner();
//    const unallocLabelList& neighbour = mesh().neighbour();
//
//    // new formulation for deltaX
//    const volVectorField& cellCenter = mesh().C();
//    const surfaceVectorField& faceCenter = mesh().Cf();
//
//    // Get the face area vector
//    const surfaceVectorField& Sf = mesh_.Sf();      //face area vector
//    const surfaceScalarField& magSf = mesh_.magSf();//face area
//
//    volScalarField deltaX
//    (
//        IOobject
//        (
//            "deltaXEdge",
//            mesh().time().timeName(),
//            mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("deltaXEdge", dimLength, HUGE),
//        zeroGradientFvPatchScalarField::typeName
//    );
//
//    forAll(owner, faceI)
//    {
//        label own = owner[faceI];
//        label nei = neighbour[faceI];
//
//        vector normal = Sf[faceI]/magSf[faceI]; //face normal vector
//
//        deltaX[own] =           //why for each owner + neighbor cell???
//            min(deltaX[own], mag((faceCenter[faceI]-cellCenter[own]) & normal));
//
//        deltaX[nei] =
//            min(deltaX[nei], mag((faceCenter[faceI]-cellCenter[nei]) & normal));
//    }
//
//    forAll(deltaX.boundaryField(), patchI)
//    {
//
//        const fvPatchScalarField& pp = deltaX.boundaryField()[patchI];
//        const vectorField pFaceCenter = pp.patch().Cf();
//
//        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
//        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];
//
//        const unallocLabelList& faceCells =  pp.patch().faceCells();
//
//
//        forAll(pp, facei)
//        {
//            vector normal = pSf[facei]/pMagSf[facei];
//
//            label own = faceCells[facei];
//
//            deltaX[own] = min
//            (
//                deltaX[own],
//                mag((pFaceCenter[facei]-cellCenter[own]) & normal)
//            );
//        }
//    }
//
//    deltaX.correctBoundaryConditions();
////     deltaX.write();
//    deltaX_ = deltaX;
//    deltaX_.correctBoundaryConditions();
//
//};


void localTimeStep::update(scalar maxCo, Switch adjustTimeStep)
{
//    if (mesh().moving())
//    {
//        updateDeltaX();
//    }

    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();
//    // Get the face area vector
//    const surfaceVectorField& Sf = mesh_.Sf();
//    // Get the face areas
//    const surfaceScalarField& magSf = mesh_.magSf();


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



    surfaceScalarField deltaCoeffs = mesh_.surfaceInterpolation::deltaCoeffs();

    //inviscid deltaT on the faces
 //   surfaceScalarField deltaT_invf = 1/((mag(Uf_)+cf_)*deltaCoeffs);
    surfaceVectorField normal = Sf/magSf; //face normal vector
    surfaceScalarField deltaT_invf = 1/((mag(Uf_ & normal)+cf_)*deltaCoeffs);

    surfaceScalarField deltaT_visf =  1/(max(4/(3*rhof_),kappaf_/rhof_)
                                        *muEfff_*sqr(deltaCoeffs));

    // loop over all faces
    // find minimum deltaT at each face and move it into the cell
    forAll(Uf_, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        deltaTMinValue_inv[own] =  min(deltaTMinValue_inv[own], deltaT_invf[faceI]);
        deltaTMinValue_inv[nei] =  min(deltaTMinValue_inv[nei], deltaT_invf[faceI]);

        deltaTMinValue_vis[own] =  min(deltaTMinValue_vis[own], deltaT_visf[faceI]);
        deltaTMinValue_vis[nei] =  min(deltaTMinValue_vis[nei], deltaT_visf[faceI]);
     
    }
    
    forAll(Uf_, faceI)
    {
        label own = owner[faceI];
        //label nei = neighbour[faceI];


	lambdaMax[faceI] = max(max((Uf_[faceI] & normal[faceI])+cf_[faceI],
				(Uf_[faceI] & normal[faceI])-cf_[faceI]),(Uf_[faceI] & normal[faceI]));

	phif_[own] += (magSf[faceI]*lambdaMax[faceI]);
//	phif_[own] += (magSf[faceI]*(mag(Uf_[faceI]&normal[faceI])+cf_[faceI]));	
     
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
	    
	    
	    plambdaMax[facei]= max(max((pUf[facei] & normal)+pcf[facei],
				(pUf[facei] & normal)-pcf[facei]),(pUf[facei] & normal));

	    phif_[own] += (pMagSf[facei]*plambdaMax[facei]);					
//	    phif_[own] += (pMagSf[facei]*(mag(pUf[facei]&normal)+pcf[facei]));				

//            deltaX[own] = min
//            (
//                deltaX[own],
//                mag((pFaceCenter[facei]-cellCenter[own]) & normal)
//            );

        }
    }       

    maxDeltaTFluent_ = 2.0*cellVolume/phif_;
    

//    maxDeltaT_ = ((deltaTMinValue_inv * deltaTMinValue_inv)
//                        /(deltaTMinValue_inv + deltaTMinValue_inv));

    maxDeltaT_ = min(deltaTMinValue_inv, deltaTMinValue_vis);
      


    CoDeltaT_ = maxCo* maxDeltaTFluent_;

//    // loop over all cells
//    forAll(deltaT_inv, cellI)
//    {
////        Lambda_vis[cellI] = Lambda_vis[cellI]/cellVolume[cellI];
//
//        maxDeltaT_[cellI] = (deltaT_inv[cellI]* deltaT_vis[cellI])
//                /(deltaT_inv[cellI]+ deltaT_vis[cellI]);
//    }

    

//    volScalarField speedOfSound = thermophysicalModel_.Cp()/
//        (thermophysicalModel_.Cv()*thermophysicalModel_.psi());
//
//    dimensionedScalar minVelocity("UMin", U_.dimensions(), SMALL);
//
//    bound(speedOfSound,minVelocity);
//
//    // compute the maximum inviscid deltaT
//    volScalarField deltaTInvis
//    (
//       deltaX_/(mag(U_)+sqrt(speedOfSound))
//    );
//
//    if ( max(turbulenceModel_.muEff()).value() > SMALL )
//    {
//        // compute the maximum viscous deltaT
//        volScalarField deltaTVis
//        (
//            sqr(deltaX_)*thermophysicalModel_.rho()/turbulenceModel_.muEff()
//        );
//
//        // blend both time step values and limit them with the given
//        // maximum CFL number, has only some influence if both timesteps are in the same order
//        // otherwise it is almost the smaller one of both
//        // This timestep is always smaller than the min() of both
//        // cf. Arnone et al.
//        CoDeltaT_ = maxCo*(deltaTVis*deltaTInvis)/(deltaTVis+deltaTInvis);
//        // cf. Weiss and Smith
////         CoDeltaT_ = maxCo*min(deltaTVis,deltaTInvis);
//    }
//
////dimensionedScalar tMin("tMin", dimTime, 2.0e-7);
////    CoDeltaT_ = tMin;
//
//    else
//    {
//        CoDeltaT_ = maxCo*deltaTInvis;
//    }

    if (adjustTimeStep)
    {
        CoDeltaT_ = min(CoDeltaT_);
    }

}


void localTimeStep::updateLast(scalar maxCo, Switch adjustTimeStep)
{
    CoDeltaT_ = 1.0*min(CoDeltaT_);
}

void localTimeStep::update(scalar maxCo)
{
    update(maxCo,"no");
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
