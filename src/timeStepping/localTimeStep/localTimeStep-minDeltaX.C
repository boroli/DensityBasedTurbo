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
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("CoDeltaT", dimTime, 0.1),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaX_
    (
        IOobject
        (
            "deltaX",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaX", dimLength, SMALL),
        zeroGradientFvPatchScalarField::typeName
    ),

    deltaTVisOverDeltaTInv_
    (
        IOobject
        (
            "deltaTVisOverDeltaTInv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (sqr(deltaX_)*thermophysicalModel_.rho()/turbulenceModel_.muEff())
        /(deltaX_/(mag(U)+sqrt(thermophysicalModel_.Cp()/(thermophysicalModel_.Cv()
            *thermophysicalModel_.psi()))))
    )


{
    updateDeltaX();
};

void localTimeStep::updateDeltaX()
{
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();
/*
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

    volScalarField maxFaceArea
    (
        IOobject
        (
            "maxFaceArea",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("maxFaceArea", dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    const surfaceScalarField& faceAreas = mesh().magSf();

    // compute maximum face area for each cell from the internal faces
    forAll(owner, facei)
    {
        maxFaceArea[owner[facei]] =
            max(maxFaceArea[owner[facei]], faceAreas[facei]);

        maxFaceArea[neighbour[facei]] =
            max(maxFaceArea[neighbour[facei]], faceAreas[facei]);
    }

    // compute maximum face area for each cell from the boundary faces
    forAll(maxFaceArea.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pfaceArea =
            faceAreas.boundaryField()[patchi];

        const fvPatch& p = pfaceArea.patch();
        const unallocLabelList& faceCells = p.patch().faceCells();

        forAll(pfaceArea, patchFacei)
        {
            maxFaceArea[faceCells[patchFacei]] = max
            (
                maxFaceArea[faceCells[patchFacei]],
                pfaceArea[patchFacei]
            );
        }
    }

    // update boundary conditons for boundary faces
    maxFaceArea.correctBoundaryConditions();

    // compute characteristic length for each cell
    deltaX_ = cellVolume/maxFaceArea;
    deltaX_.correctBoundaryConditions();
*/
//     volScalarField deltaXVol(cellVolume/maxFaceArea);
//     deltaXVol.rename("deltaXVol");
//     deltaXVol.write();

    // new formulation for deltaX
    const volVectorField& cellCenter = mesh().C();
    const surfaceVectorField& faceCenter = mesh().Cf();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();      //face area vector
    const surfaceScalarField& magSf = mesh_.magSf();//face area

    volScalarField deltaX
    (
        IOobject
        (
            "deltaXEdge",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaXEdge", dimLength, HUGE),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector normal = Sf[faceI]/magSf[faceI]; //face normal vector

        deltaX[own] =           //why for each owner + neighbor cell???
            min(deltaX[own], mag((faceCenter[faceI]-cellCenter[own]) & normal));

        deltaX[nei] =
            min(deltaX[nei], mag((faceCenter[faceI]-cellCenter[nei]) & normal));
    }

    forAll(deltaX.boundaryField(), patchI)
    {

        const fvPatchScalarField& pp = deltaX.boundaryField()[patchI];
        const vectorField pFaceCenter = pp.patch().Cf();

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];

        const unallocLabelList& faceCells =  pp.patch().faceCells();

        forAll(pp, facei)
        {
            vector normal = pSf[facei]/pMagSf[facei];

            label own = faceCells[facei];

            deltaX[own] = min
            (
                deltaX[own],
                mag((pFaceCenter[facei]-cellCenter[own]) & normal)
            );
        }
    }


//     forAll(owner, faceI)
//     {
//         label own = owner[faceI];
//         label nei = neighbour[faceI];
// 
// //         scalar deltaLLeft  = mag(faceCenter[faceI] - cellCenter[own]);
// //         scalar deltaLRight = mag(faceCenter[faceI] - cellCenter[nei]);
// //         scalar deltaL = mag(fcellCenter[nei] - cellCenter[own]);
// 
//         deltaX[own] =
//             min(deltaX[own], mag(faceCenter[faceI]-cellCenter[own]));
// 
//         deltaX[nei] =
//             min(deltaX[nei], mag(faceCenter[faceI]-cellCenter[nei]));
//     }
// 
//     forAll(deltaX.boundaryField(), patchI)
//     {
//         const fvsPatchVectorField& pFaceCenter =
//             faceCenter.boundaryField()[patchI];
// 
//         const fvPatch& p = pFaceCenter.patch();
//         const unallocLabelList& faceCells = p.patch().faceCells();
// 
//         forAll(pFaceCenter, patchFaceI)
//         {
//             label cellI = faceCells[patchFaceI];
//             deltaX[cellI] = min
//             (
//                 deltaX[cellI],
//                 mag(pFaceCenter[patchFaceI]-cellCenter[cellI])
//             );
//         }
//     }

    deltaX.correctBoundaryConditions();
//     deltaX.write();
    deltaX_ = deltaX;
    deltaX_.correctBoundaryConditions();

};


void localTimeStep::update(scalar maxCo, Switch adjustTimeStep)
{
    if (mesh().moving())
    {
        updateDeltaX();
    }

    volScalarField speedOfSound = thermophysicalModel_.Cp()/
        (thermophysicalModel_.Cv()*thermophysicalModel_.psi());

    dimensionedScalar minVelocity("UMin", U_.dimensions(), SMALL);

    bound(speedOfSound,minVelocity);

    // compute the maximum inviscid deltaT
    volScalarField deltaTInvis
    (
       deltaX_/(mag(U_)+sqrt(speedOfSound))
    );

    if ( max(turbulenceModel_.muEff()).value() > SMALL )
    {
        // compute the maximum viscous deltaT
        volScalarField deltaTVis
        (
            sqr(deltaX_)*thermophysicalModel_.rho()/turbulenceModel_.muEff()
        );

        // blend both time step values and limit them with the given
        // maximum CFL number, has only some influence if both timesteps are in the same order
        // otherwise it is almost the smaller one of both
        // This timestep is always smaller than the min() of both
        // cf. Arnone et al.
         CoDeltaT_ = maxCo*(deltaTVis*deltaTInvis)/(deltaTVis+deltaTInvis);
        // cf. Weiss and Smith
//         CoDeltaT_ = maxCo*min(deltaTVis,deltaTInvis);

    deltaTVisOverDeltaTInv_ = deltaTVis/deltaTInvis;
    }
    
//dimensionedScalar tMin("tMin", dimTime, 2.0e-7);
//    CoDeltaT_ = tMin;

    else
    {
        CoDeltaT_ = maxCo*deltaTInvis;
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
