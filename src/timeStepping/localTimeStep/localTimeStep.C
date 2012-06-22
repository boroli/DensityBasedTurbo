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
    deltaS_
    (
        IOobject
        (
            "deltaS",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("deltaS", dimArea, vector::one),
        zeroGradientFvPatchVectorField::typeName
    ),
    cellVolume_
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
        dimensionedScalar("cellVolume", dimVolume, 0.1),
        zeroGradientFvPatchScalarField::typeName
    )
{
    updateDeltaX();
};

void localTimeStep::updateDeltaX()
{
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    // new formulation for deltaX according to Blazek

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();

    // Reset values
    deltaS_ = dimensionedVector("deltaS", dimArea, vector::zero);
    cellVolume_.internalField() = mesh().V();
    cellVolume_.correctBoundaryConditions();

    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        deltaS_[own] += 0.5*cmptMag(Sf[faceI]);
        deltaS_[nei] += 0.5*cmptMag(Sf[faceI]);
    }

    forAll(deltaS_.boundaryField(), patchI)
    {
        const fvPatchVectorField& pp = deltaS_.boundaryField()[patchI];

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];

        const unallocLabelList& faceCells =  pp.patch().faceCells();

        forAll(pp, faceI)
        {
            label own = faceCells[faceI];

            deltaS_[own] += 0.5*cmptMag(pSf[faceI]);
        }
    }

    deltaS_.correctBoundaryConditions();
};

void localTimeStep::update(scalar maxCo, Switch adjustTimeStep)
{
    if (mesh().moving())
    {
        updateDeltaX();
    }
    
    // square of speed of sound for ideal gases
    volScalarField speedOfSound = thermophysicalModel_.Cp()/
        (thermophysicalModel_.Cv()*thermophysicalModel_.psi());
    dimensionedScalar minVelocity("UMin", U_.dimensions(), SMALL);
    bound(speedOfSound,minVelocity);

    // constant for viscous time step is hard coded for the moment
    scalar constantC = 4.0;

    if ( max(turbulenceModel_.muEff()).value() > SMALL )
    {
        CoDeltaT_
            = maxCo*cellVolume_/(((U_+sqrt(speedOfSound)*vector::one) & deltaS_) 
            + constantC
            *max(4./3.,thermophysicalModel_.Cp()/thermophysicalModel_.Cv())
            *turbulenceModel_.alphaEff()
            *magSqr(deltaS_)/(thermophysicalModel_.rho()*cellVolume_) );
    }
    else
    {
        CoDeltaT_ 
            = maxCo*cellVolume_/((U_+sqrt(speedOfSound)*vector::one) & deltaS_);
    }

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
