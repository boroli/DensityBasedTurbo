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
    godunovFlux

Description
    generic Godunov flux class

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "godunovFlux.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Flux> // For the class
template<class phiType, class gradPhiType, class Limiter> // For the function
void Foam::godunovFlux<Flux>::updateLimiter
(
    const GeometricField<phiType, fvPatchField, volMesh>& phi,
    const GeometricField<gradPhiType, fvPatchField, volMesh>& gradPhi,
    GeometricField<phiType, fvPatchField, volMesh>& phiLimiter,
    word oneDLimiterName
)
{
    // Reset limiter field
    phiLimiter = pTraits<phiType>::one;
    phiLimiter.correctBoundaryConditions();

    if (multidimLimiterSwitch_)
    {
        // 2nd order correction
        IStringStream blendingFactor(name(epsilon_));

        Limiter phiLimiterFunction(blendingFactor);

        GeometricField<phiType, fvPatchField, volMesh> phiMinValue("phiMin",phi);
        GeometricField<phiType, fvPatchField, volMesh> phiMaxValue("phiMax",phi);

        // for BJ/VK find min/max for each phi
        // internal face
        forAll(owner_, faceI)
        {
            label own = owner_[faceI];
            label nei = neighbour_[faceI];

            // min value owner
            phiMinValue[own] =
                min(phiMinValue[own], phi[nei]);
            // max value owner
            phiMaxValue[own] =
                max(phiMaxValue[own], phi[nei]);

            // min value neighbour
            phiMinValue[nei] =
                min(phiMinValue[nei], phi[own]);
            // max value neighbour
            phiMaxValue[nei] =
                max(phiMaxValue[nei], phi[own]);
        }

        // Update coupled boundary min/max values of primitive variables
        forAll(phi.boundaryField(), patchi)
        {
            const fvPatchField<phiType>& pphi = phi.boundaryField()[patchi];
            const unallocLabelList& faceCells = pphi.patch().faceCells();

            if (pphi.coupled())
            {
                // primitive variables
                const Field<phiType> pphiRight = pphi.patchNeighbourField();

                forAll(pphi, faceI)
                {
                    label own = faceCells[faceI];

                    // min values at coupled boundary faces
                    phiMinValue[own] =
                        min(phiMinValue[own], pphiRight[faceI]);

                    // max values at coupled boundary faces
                    phiMaxValue[own] =
                        max(phiMaxValue[own], pphiRight[faceI]);
                }
            }
        }

        // o.b. An update of the boundary conditions for
        // the min/max values seems to be needed for coupled patches
        phiMinValue.correctBoundaryConditions();
        phiMaxValue.correctBoundaryConditions();

        // compute for each cell a limiter
        // Loop over all faces with different deltaR vector
        forAll(owner_, faceI)
        {
            label own = owner_[faceI];
            label nei = neighbour_[faceI];

            // find minimal limiter value in each cell
            phiLimiter[own] =
            min
            (
                phiLimiter[own], phiLimiterFunction.limiter
                (
                    cellVolume_[own],
                    1.0, // flux dummy
                    phiMaxValue[own]-phi[own],
                    phiMinValue[own]-phi[own],
                    gradPhi[own],
                    gradPhi[own],
                    faceCenter_[faceI]-cellCenter_[own]
                )
            );

            phiLimiter[nei] =
            min
            (
                phiLimiter[nei],
                phiLimiterFunction.limiter
                (
                    cellVolume_[nei],
                    1.0, // flux dummy
                    phiMaxValue[nei]-phi[nei],
                    phiMinValue[nei]-phi[nei],
                    gradPhi[nei],
                    gradPhi[nei],
                    faceCenter_[faceI]-cellCenter_[nei]
                )
            );
        }

        // Update coupled boundary limiters
        forAll(phi.boundaryField(), patchi)
        {
            const fvPatchField<phiType>& pphi = phi.boundaryField()[patchi];
            const unallocLabelList& faceCells = pphi.patch().faceCells();

            if (pphi.coupled())
            {
                // cell and face centers
                const vectorField delta =  pphi.patch().delta();

                forAll(pphi, faceI)
                {
                    label own = faceCells[faceI];

                    phiLimiter[own] = 
                    min
                    (
                        phiLimiter[own],
                        phiLimiterFunction.limiter
                        (
                            cellVolume_[own],
                            1.0, // flux dummy
                            phiMaxValue[own]-phi[own],
                            phiMinValue[own]-phi[own],
                            gradPhi[own],
                            gradPhi[own],
                            delta[faceI]
                        )
                    );
                }
            }
//             else
//             {
//                 forAll(pphi, faceI)
//                 {
//                     label own = faceCells[faceI];
//                     // Calculate minimal limiters
//                     phiLimiter[own] =
//                         min(phiLimiter[own], pTraits<phiType>::one);
//                 }
//             }
        }
    }
    else
    {
        IStringStream phiSchemeData(oneDLimiterName);

        // compute face based limiters
        surfaceScalarField phiLimiterSurface =
            limitedSurfaceInterpolationScheme<phiType>::New
            (
                mesh_,
                rhoFlux_,
                phiSchemeData
            )().limiter(phi);

        // use smallest face based limiter for volume
        forAll(owner_, faceI)
        {
            label own = owner_[faceI];
            label nei = neighbour_[faceI];

            // find minimal limiter value in each cell
            phiLimiter[own] =
                min(phiLimiter[own], phiLimiterSurface[faceI]*pTraits<phiType>::one);

            phiLimiter[nei] =
                min(phiLimiter[nei], phiLimiterSurface[faceI]*pTraits<phiType>::one);
        }

        // Update coupled boundary limiters
        forAll(phi.boundaryField(), patchi)
        {
            const fvPatchField<phiType>& pphi = phi.boundaryField()[patchi];
            const unallocLabelList& faceCells = pphi.patch().faceCells();

            if (pphi.coupled())
            {
                forAll(pphi, faceI)
                {
                    label own = faceCells[faceI];

                    phiLimiter[own] =
                        min(phiLimiter[own], phiLimiterSurface[faceI]*pTraits<phiType>::one);
                }
            }
        }
    }

    // o.b. An update of the boundary conditions for
    // the limiter values seems to be needed for coupled patches
    phiLimiter.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux>
Foam::godunovFlux<Flux>::godunovFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& rho,
    const basicThermo& thermophysicalModel,
    const compressible::turbulenceModel& turbulenceModel
)
:
    mesh_(p.mesh()),
    p_(p),
    U_(U),
    rho_(rho),
    thermophysicalModel_(thermophysicalModel),
    turbulenceModel_(turbulenceModel),
    owner_(mesh_.owner()),
    neighbour_(mesh_.neighbour()),
    cellCenter_(mesh_.C()),
    faceCenter_(mesh_.Cf()),
    k_(turbulenceModel_.k()),
    kappa_(thermophysicalModel_.Cp()/thermophysicalModel_.Cv()),
    cellVolume_
    (
        IOobject
        (
            "cellVolume",
            mesh_.time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimVolume,
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoFlux_
    (
        IOobject
        (
            "rhoFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // fluxes in mass conservation equation: \varrho \vec{u}
        // only initialisation!
        (linearInterpolate(rho_*U_) & mesh_.Sf())
    ),
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // fluxes in momentum equation: \varrho \vec{u} \vec{u}
        // only initialisation!
        rhoFlux_*linearInterpolate(U_)
    ),
    rhoEFlux_
    (
        IOobject
        (
            "rhoEFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // fluxes in total energy equation: \varrho E \vec{u}
        // only initialisation!
        rhoFlux_*linearInterpolate
        (
            thermophysicalModel_.h() + 0.5*magSqr(U_) + k_ - (p_/rho_)
        )
    ),
    dotX_
    (
        IOobject
        (
            "dotX",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimVelocity, vector::zero)
    ),
    gradp_(fvc::grad(p_,"grad(pSlope)")),
    gradU_(fvc::grad(U_,"grad(USlope)")),
    gradrho_(fvc::grad(rho_,"grad(rhoSlope)")),
    gradk_(fvc::grad(k_,"grad(TKE)")),
    pLimiter_
    (
        IOobject
        (
            "pLimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("pLimiter", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    ULimiter_
    (
        IOobject
        (
            "ULimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("ULimiter", dimless, vector::one),
        zeroGradientFvPatchVectorField::typeName
    ),
    rhoLimiter_
    (
        IOobject
        (
            "rhoLimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("rhoLimiter", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    kLimiter_
    (
        IOobject
        (
            "kLimiter",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("kLimiter", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    epsilon_(5),
    Konstant_(0.05),
    limiterName_("vanAlbadaSlope"),
    multidimLimiterSwitch_(false)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux>
void Foam::godunovFlux<Flux>::update(Switch secondOrder)
{
    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    // update cell volume
    cellVolume_.internalField() = mesh_.V();
    cellVolume_.correctBoundaryConditions();

    // need a copies here, because the return value is a <tmp>
    k_ = turbulenceModel_.k();
    // WARNING: only valid for ideal gases!
    kappa_ = thermophysicalModel_.Cp()/thermophysicalModel_.Cv();

    // read riemann solver coeffs
    if(mesh_.solutionDict().found("Riemann"))
    {
        // Read Riemann coeffs dict
        dictionary riemann = mesh_.solutionDict().subDict("Riemann");
        if (riemann.found("multidimLimiter"))
        {
            multidimLimiterSwitch_ = Switch(riemann.lookup("multidimLimiter"));
        }
        if (riemann.found("limiterName"))
        {
            limiterName_ = word(riemann.lookup("limiterName"));
        }
        epsilon_ = riemann.lookupOrDefault("epsilon",epsilon_);
        Konstant_ = riemann.lookupOrDefault("RiemannSolverKonstant",Konstant_);
    }

    // 2nd order correction
    if (secondOrder)
    {
        // Update the primitive gradients
        gradp_ = fvc::grad(p_,"grad(pSlope)");
        gradU_ = fvc::grad(U_,"grad(USlope)");
        gradrho_ = fvc::grad(rho_,"grad(rhoSlope)");
        if (max(k_).value() > 0.0)
        {
            gradk_ = fvc::grad(k_,"grad(TKE)");
        }

        // MultiDimensional Limiters
//         updateLimiter<scalar,vector,BarthJespersenSlopeMultiLimiter>(p_,gradp_,pLimiter_);
// //
        updateLimiter<scalar,vector,VenkatakrishnanSlopeMultiLimiter>(p_,gradp_,pLimiter_,limiterName_);
        updateLimiter<vector,tensor,VenkatakrishnanSlopeMultiLimiter>(U_,gradU_,ULimiter_,limiterName_+"V");
        updateLimiter<scalar,vector,VenkatakrishnanSlopeMultiLimiter>(rho_,gradrho_,rhoLimiter_,limiterName_);
        updateLimiter<scalar,vector,VenkatakrishnanSlopeMultiLimiter>(k_,gradk_,kLimiter_,limiterName_);
    }

//     boundMinMax(pLimiter_,dimensionedScalar(0.0),dimensionedScalar(1.0));
//     boundMinMax(ULimiter_[0],dimensionedScalar(0.0),dimensionedScalar(1.0));
//     boundMinMax(rhoLimiter,dimensionedScalar(0.0),dimensionedScalar(1.0));

//     Info << "max(pLimiter_) "  << max(pLimiter_.internalField())
//          << " min(pLimiter_) " << min(pLimiter_.internalField())  << endl;
//     Info << "max(ULimiter_) "  << max(ULimiter_.internalField())
//          << " min(ULimiter_) " << min(ULimiter_.internalField())  << endl;
//     Info << "max(rhoLimiter_) "  << max(rhoLimiter_.internalField())
//          << " min(rhoLimiter_) " << min(rhoLimiter_.internalField())  << endl;
//     Info << "max(kLimiter_) "  << max(kLimiter_.internalField())
//          << " min(kLimiter_) " << min(kLimiter_.internalField())  << endl;

    // Calculate fluxes at internal faces
    forAll(owner_, faceI)
    {
        label own = owner_[faceI];
        label nei = neighbour_[faceI];

        vector deltaRLeft  = faceCenter_[faceI] - cellCenter_[own];
        vector deltaRRight = faceCenter_[faceI] - cellCenter_[nei];

        // calculate fluxes with reconstructed primitive variables at faces
        // TODO: thermophysical variables are not reconstructed at faces!!!
        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            p_[own] + secondOrder*pLimiter_[own]*(deltaRLeft  & gradp_[own]), // reconstructed left p
            p_[nei] + secondOrder*pLimiter_[nei]*(deltaRRight & gradp_[nei]), // reconstructed right p
            U_[own] + secondOrder*cmptMultiply(ULimiter_[own],(deltaRLeft  & gradU_[own])), // reconstructed left U
            U_[nei] + secondOrder*cmptMultiply(ULimiter_[nei],(deltaRRight & gradU_[nei])), // reconstructed right U
//
//          using minimum component scalar limiter
//             U_[own] + secondOrder*cmptMin(ULimiter_[own])*(deltaRLeft & gradU_[own]),  // reconstructed left U
//             U_[nei] + secondOrder*cmptMin(ULimiter_[nei])*(deltaRRight & gradU_[nei]), // reconstructed right U
//
            rho_[own] + secondOrder*rhoLimiter_[own]*(deltaRLeft  & gradrho_[own]), // reconstructed left rho
            rho_[nei] + secondOrder*rhoLimiter_[nei]*(deltaRRight & gradrho_[nei]), // reconstructed right rho
            k_[own] + secondOrder*kLimiter_[own]*(deltaRLeft  & gradk_[own]),
            k_[nei] + secondOrder*kLimiter_[nei]*(deltaRRight & gradk_[nei]),
// //
            kappa_[own],       // left kappa
            kappa_[nei],       // right kappa
            Sf[faceI],      // face vector
            magSf[faceI],   // face area
            dotX_[faceI],    // face velocity
            Konstant_
        );
    }

    // Update boundary field and values
    forAll(p_.boundaryField(), patchi)
    {
        fvsPatchScalarField& pRhoFlux = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryField()[patchi];

        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const fvPatchVectorField& pU = U_.boundaryField()[patchi];
        const fvPatchScalarField& prho = rho_.boundaryField()[patchi];
        const fvPatchScalarField& pk = k_.boundaryField()[patchi];

        const fvPatchVectorField& pGradp = gradp_.boundaryField()[patchi];
        const fvPatchTensorField& pGradU = gradU_.boundaryField()[patchi];
        const fvPatchVectorField& pGradrho = gradrho_.boundaryField()[patchi];
        const fvPatchVectorField& pGradk = gradk_.boundaryField()[patchi];

        const fvPatchScalarField& pkappa = kappa_.boundaryField()[patchi];

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const fvsPatchVectorField& pDotX = dotX_.boundaryField()[patchi];

        const fvPatchVectorField& pCellCenter = cellCenter_.boundaryField()[patchi];
//         const fvsPatchVectorField& pFaceCenter = faceCenter.boundaryField()[patchi];

        const fvPatchScalarField& ppLimiter = pLimiter_.boundaryField()[patchi];
        const fvPatchVectorField& pULimiter = ULimiter_.boundaryField()[patchi];
        const fvPatchScalarField& prhoLimiter = rhoLimiter_.boundaryField()[patchi];
        const fvPatchScalarField& pkLimiter = kLimiter_.boundaryField()[patchi];

        // special treatment at coupled boundaries
        if (pp.coupled())
        {
            // primitive variables
            const scalarField ppLeft  = pp.patchInternalField();
            const scalarField ppRight = pp.patchNeighbourField();

            const vectorField pULeft  = pU.patchInternalField();
            const vectorField pURight = pU.patchNeighbourField();

            const scalarField prhoLeft  = prho.patchInternalField();
            const scalarField prhoRight = prho.patchNeighbourField();

            const scalarField pkLeft  = pk.patchInternalField();
            const scalarField pkRight = pk.patchNeighbourField();

            const scalarField pkappaLeft  = pkappa.patchInternalField();
            const scalarField pkappaRight = pkappa.patchNeighbourField();

            // cell gradients
            const vectorField pGradpLeft  = pGradp.patchInternalField();
            const vectorField pGradpRight = pGradp.patchNeighbourField();

            const tensorField pGradULeft  = pGradU.patchInternalField();
            const tensorField pGradURight = pGradU.patchNeighbourField();

            const vectorField pGradrhoLeft  = pGradrho.patchInternalField();
            const vectorField pGradrhoRight = pGradrho.patchNeighbourField();

            const vectorField pGradkLeft  = pGradk.patchInternalField();
            const vectorField pGradkRight = pGradk.patchNeighbourField();

            // cell limiters
            scalarField ppLimiterLeft  = ppLimiter.patchInternalField();
            scalarField ppLimiterRight = ppLimiter.patchNeighbourField();

            vectorField pULimiterLeft  = pULimiter.patchInternalField();
            vectorField pULimiterRight = pULimiter.patchNeighbourField();

            scalarField prhoLimiterLeft  = prhoLimiter.patchInternalField();
            scalarField prhoLimiterRight = prhoLimiter.patchNeighbourField();

            scalarField pkLimiterLeft  = pkLimiter.patchInternalField();
            scalarField pkLimiterRight = pkLimiter.patchNeighbourField();

            // cell and face centers
            const vectorField faceCenter =  pp.patch().Cf();
//             const vectorField deltaRLeftField = pp.patch().delta();
            const vectorField pCellCenterLeft  =  pCellCenter.patchInternalField();
            const vectorField pCellCenterRight =  pCellCenter.patchNeighbourField();

            forAll(pp, faceI)
            {
                vector deltaRLeft  = faceCenter[faceI] - pCellCenterLeft[faceI];
                vector deltaRRight = faceCenter[faceI] - pCellCenterRight[faceI];

                // bound the limiters between 0 and 1 in order to prevent problems due to interpolations
                ppLimiterLeft[faceI]   = max(min(ppLimiterLeft[faceI],1.0),0.0);
                ppLimiterRight[faceI]  = max(min(ppLimiterRight[faceI],1.0),0.0);

                pULimiterLeft[faceI]  = max(min(pULimiterLeft[faceI],vector::one),vector::zero);
                pULimiterRight[faceI] = max(min(pULimiterRight[faceI],vector::one),vector::zero);

                prhoLimiterLeft[faceI]   = max(min(prhoLimiterLeft[faceI],1.0),0.0);
                prhoLimiterRight[faceI]  = max(min(prhoLimiterRight[faceI],1.0),0.0);

                pkLimiterLeft[faceI]   = max(min(pkLimiterLeft[faceI],1.0),0.0);
                pkLimiterRight[faceI]  = max(min(pkLimiterRight[faceI],1.0),0.0);

                // Calculate fluxes at coupled boundary faces
                Flux::evaluateFlux
                (
                    pRhoFlux[faceI],
                    pRhoUFlux[faceI],
                    pRhoEFlux[faceI],
                    ppLeft[faceI]  + secondOrder*ppLimiterLeft[faceI] *(deltaRLeft  & pGradpLeft[faceI]),                // face p
                    ppRight[faceI] + secondOrder*ppLimiterRight[faceI]*(deltaRRight & pGradpRight[faceI]),               // face p
                    pULeft[faceI]  + secondOrder*cmptMultiply(pULimiterLeft[faceI] ,(deltaRLeft  & pGradULeft[faceI])),  // face U
                    pURight[faceI] + secondOrder*cmptMultiply(pULimiterRight[faceI],(deltaRRight & pGradURight[faceI])), // face U
//
//                     pULeft[faceI]  + secondOrder*cmptMin(pULimiterLeft[faceI])*(deltaRLeft & pGradULeft[faceI]),      // face U
//                     pURight[faceI] + secondOrder*cmptMin(pULimiterRight[faceI])*(deltaRRight & pGradURight[faceI]),   // face U
//
                    prhoLeft[faceI]  + secondOrder*prhoLimiterLeft[faceI] *(deltaRLeft  & pGradrhoLeft[faceI]),          // face rho
                    prhoRight[faceI] + secondOrder*prhoLimiterRight[faceI]*(deltaRRight & pGradrhoRight[faceI]),         // face rho
                    pkLeft[faceI]  + secondOrder*pkLimiterLeft[faceI] *(deltaRLeft  & pGradkLeft[faceI]),                // face k
                    pkRight[faceI] + secondOrder*pkLimiterRight[faceI]*(deltaRRight & pGradkRight[faceI]),               // face k
// //
                    pkappaLeft[faceI],  // face kappa
                    pkappaRight[faceI], // face kappa
                    pSf[faceI],         // face vector
                    pMagSf[faceI],      // face area
                    pDotX[faceI],       // face velocity
                    Konstant_
                );
            }
        }
        else
        {
            forAll(pp, faceI)
            {
                // Calculate fluxes at boundary faces
                Flux::evaluateFlux
                (
                    pRhoFlux[faceI],
                    pRhoUFlux[faceI],
                    pRhoEFlux[faceI],
                    pp[faceI],        // face p
                    pp[faceI],        // face p
                    pU[faceI],        // face U
                    pU[faceI],        // face U
                    prho[faceI],      // face rho
                    prho[faceI],      // face rho
                    pk[faceI],        // face k
                    pk[faceI],        // face k
// //
                    pkappa[faceI],    // face kappa
                    pkappa[faceI],    // face kappa
                    pSf[faceI],       // face vector
                    pMagSf[faceI],    // face area
                    pDotX[faceI],     // face velocity
                    Konstant_
                );
            }
        }
    }
}

// ************************************************************************* //
