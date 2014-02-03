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
    roeMultiDensityALEFlux

Description

Author
    Hrvoje Jasak All rights reserved.
    Aleksandar Jemcov  All rights reserved.
    Oliver Borm  All rights reserved.
    Sebastian Saegeler All rights reserved.

\*---------------------------------------------------------------------------*/

#include "roeMultiDensityALEFlux.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::roeMultiDensityALEFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,
// //
//     const scalar TLeft,
//     const scalar TRight,
// //
    const scalar rhoLeft,
    const scalar rhoRight,
// //
    const scalar visVelLeft,
    const scalar visVelRight,
// //
     const scalar RLeft,
     const scalar RRight,
//     const scalar CvLeft,
//     const scalar CvRight,
// //
    const scalar kappaLeft,
    const scalar kappaRight,
// //
    const vector Sf,
    const scalar magSf,
    const vector dotX,
    const scalar Konstant       
) const
{
//    scalar R = 287.0;    
//    scalar Cp = 1007.0;

    // Decode left and right total energy:
    const scalar eLeft =
        pLeft /(rhoLeft* (kappaLeft -1.0))+ 0.5*magSqr(ULeft);// +kLeft;
    const scalar eRight =
        pRight/(rhoRight*(kappaRight-1.0))+0.5*magSqr(URight);// +kRight;

    // normal vector
    vector normalVector = Sf/magSf;
    scalar lengthNV = mag(normalVector);

    // Compute left and right contravariant velocities:
    const scalar contrVLeft  = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);

    // Compute left and right total enthalpies:
    const scalar hLeft  = eLeft  + pLeft/rhoLeft;
    const scalar hRight = eRight + pRight/rhoRight;
 
     // Compute left and right temperatures:
//    const scalar TLeft  = pLeft/(rhoLeft*R);
//    const scalar TRight = pRight/(rhoRight*R);   
    const scalar TLeft  = pLeft/(rhoLeft*RLeft);
    const scalar TRight = pRight/(rhoRight*RRight);
    
    // Step 2: compute Roe averged quantities for face:
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight,0.0));

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft,0.0));
    const scalar rhoRightSqrt = sqrt(max(rhoRight,0.0));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1.0 - wLeft;

    const vector UTilde = ULeft*wLeft + URight*wRight;
    const scalar hTilde = hLeft*wLeft + hRight*wRight;
    const scalar TTilde = TLeft*wLeft + TRight*wRight;
    const scalar qTilde = mag(UTilde);
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;
    const scalar visVelTilde = visVelLeft*wLeft + visVelRight*wRight;
    const scalar RTilde = RLeft*wLeft + RRight*wRight;
 
    // Speed of sound
    const scalar cTilde = sqrt(max((kappaTilde - 1.0)*(hTilde - 0.5*qTildeSquare),0.0));
		
	
    // Roe averaged contravariant velocity
    const scalar contrVTilde = (UTilde & normalVector);
    
    // compute the relative velocity			
    scalar Ur = min(max(max(0.1*cTilde,qTilde),visVelTilde),cTilde);
//    scalar Ur = min(max(max(0.0001*cTilde,qTilde),visVelTilde),cTilde);
//    scalar Ur = min(max(0.0001*cTilde,qTilde),cTilde);	    				
	
    // compute the parameters for the eigenvalues
    const scalar DrhoDpTilde = kappaTilde/max(sqr(cTilde),SMALL);
    const scalar DrhoDTTilde = -rhoTilde/max(TTilde,SMALL);
    
    const scalar beta = 1.0/sqr(cTilde);  // only valid for ideal gas
//    const scalar beta  = DrhoDpTilde + DrhoDTTilde/max(rhoTilde*Cp,SMALL);
    const scalar alpha = (1.0-beta*sqr(Ur))/2.0;  // valid for ideal gas
    const scalar cDash = sqrt(sqr(alpha*contrVTilde)+sqr(Ur));
    const scalar uDash = contrVTilde*(1.0-alpha);

    // Step 3: compute some differences:
    const scalar deltaP   = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU   = URight - ULeft;
    const scalar deltaContrV = deltaU & normalVector;
    const vector deltaRhoU   = rhoRight*URight - rhoLeft*ULeft;


    // Step 4: compute eigenvalues
    // derived from algebra by hand, only for Euler equation useful
    
    const scalar Urot = dotX & normalVector;   
    scalar lambda1 = mag(uDash - cDash - Urot);
    scalar lambda2 = mag(contrVTilde - Urot);
    scalar lambda3 = mag(uDash + cDash - Urot);
    
    // Step 5: check for Harten entropy correction
    //adjustable parameter
    const scalar eps = Konstant*max(lambda1,lambda3);
    
    //const scalar eps2 = 0.2*max(lambda1,lambda3);

    if(lambda1 < eps)
    {
        lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
    }

    if(lambda2 < eps)
    {
        lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
    }

    if(lambda3 < eps)
    {
        lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
    }    

    // Step 6: compute some parameters:
//    const scalar MStar = (mag(uDash + cDash) - mag(uDash - cDash))/max(2.0*cDash,SMALL);
//    const scalar cStar = (mag(uDash + cDash) + mag(uDash - cDash))/2.0;
//    const scalar dp = MStar*deltaP + (cStar-mag(contrVTilde)+alpha*contrVTilde*MStar)*rhoTilde*deltaContrV;
//    const scalar dU = MStar*deltaContrV + (cStar -(1.0-2.0*alpha)*mag(contrVTilde)-alpha*contrVTilde*MStar)*deltaP/max(rhoTilde*sqr(Ur),SMALL);

// not sure about Harten's Entorpy-Fix here    
    const scalar MStar = (mag(lambda3) - mag(lambda1))/max(2.0*cDash,SMALL);
    const scalar cStar = (mag(lambda3) + mag(lambda1))/2.0;
    const scalar dp = MStar*deltaP + (cStar-mag(lambda2)+alpha*contrVTilde*MStar)*rhoTilde*deltaContrV;
    const scalar dU = MStar*deltaContrV + (cStar -(1.0-2.0*alpha)*mag(lambda2)-alpha*contrVTilde*MStar)*deltaP/max(rhoTilde*sqr(Ur),SMALL);
    

    // Step 7: compute the parts for the dissipation term

    // rho row:
//    const scalar D1rho = mag(contrVTilde) * deltaRho;
    const scalar D1rho = mag(lambda2) * deltaRho;    
    const scalar D2rho = dU * rhoTilde;
    const scalar D3rho = pTraits<scalar>::zero;

    // first U column
//    const vector D1U = mag(contrVTilde)*deltaRhoU;
    const vector D1U = mag(lambda2)*deltaRhoU;    

    // second U column
    const vector D2U = dU * rhoTilde * UTilde;

    // third U column
    const vector D3U = dp * normalVector;

    // E row
//    const scalar D1e = mag(contrVTilde)*(rhoRight*eRight - rhoLeft*eLeft);
    const scalar D1e = mag(lambda2)*(rhoRight*eRight - rhoLeft*eLeft);    
    const scalar D2e = dU * rhoTilde * hTilde;
    const scalar D3e = dp * contrVTilde;

    // Step 8: compute left and right fluxes

    // Left flux-vector
    const scalar fluxLeft11  = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15  = hLeft*fluxLeft11;

    // Right flux-vector
    const scalar fluxRight11  = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15  = hRight*fluxRight11;

    // Step 10: compute face flux 5-vector
    const scalar flux1 =
        0.5*(fluxLeft11  + fluxRight11  - (D1rho + D2rho + D3rho));
    const vector flux24 =
        0.5*(fluxLeft124 + fluxRight124 - (D1U + D2U + D3U));
    const scalar flux5 =
        0.5*(fluxLeft15  + fluxRight15  - (D1e + D2e + D3e));

    // Compute private data
    rhoFlux  = flux1 *magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5 *magSf;
}

// ************************************************************************* //
