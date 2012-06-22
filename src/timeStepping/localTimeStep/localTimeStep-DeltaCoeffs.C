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
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("CoDeltaT", dimTime, 0.1),
        zeroGradientFvPatchScalarField::typeName
    ),
    maxDeltaT_
    (
        IOobject
        (
            "maxDeltaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("maxDeltaT", dimTime,1.0),
        zeroGradientFvPatchVectorField::typeName
    )    
{};


void localTimeStep::update(scalar maxCo, Switch adjustTimeStep)
{

    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();


    // Velocity on the cell faces
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
    // Speed of sound on the cell faces    
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
   // Density on the cell faces 
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
    // Effective viscosity on the cell faces     
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
    // Specific heat ratio on the cell faces     
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
    // Minimum deltaT in each cell    
    volScalarField deltaTMinValue_
    (
        IOobject
        (
            "deltaTMinValue",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaTMinValue", dimTime,1.0),
        zeroGradientFvPatchVectorField::typeName
    );
	    
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();	    	   

    // Get deltaCoeffitients
    surfaceScalarField deltaCoeffs = mesh_.surfaceInterpolation::deltaCoeffs();

    //deltaT on the faces
//    surfaceScalarField deltaT_invf = 1/((mag(Uf_)+cf_)*deltaCoeffs);
    surfaceVectorField normal = Sf/magSf; //face normal vector
    surfaceScalarField deltaTf_ = 1/((mag(Uf_ & normal)+cf_)*deltaCoeffs);


    // loop over all faces
    // find minimum deltaT at each face and move it into the cell
    forAll(Uf_, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        deltaTMinValue_[own] =  min(deltaTMinValue_[own], deltaTf_[faceI]);
        deltaTMinValue_[nei] =  min(deltaTMinValue_[nei], deltaTf_[faceI]);
     
    }

    maxDeltaT_ = deltaTMinValue_;
      

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
