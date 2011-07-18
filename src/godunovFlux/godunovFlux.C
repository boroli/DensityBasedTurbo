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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux>
Foam::godunovFlux<Flux>::godunovFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const basicThermo& thermophysicalModel,
    const compressible::turbulenceModel& turbulenceModel
)
:
    mesh_(p.mesh()),
    p_(p),
    U_(U),
    T_(T),
    thermophysicalModel_(thermophysicalModel),
    turbulenceModel_(turbulenceModel),
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
        (linearInterpolate(p/(( thermophysicalModel_.Cp()
        - thermophysicalModel_.Cv() )*T)*U) & mesh_.Sf())
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
        rhoFlux_*linearInterpolate(U)
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
        rhoFlux_*linearInterpolate(thermophysicalModel_.Cv()*T+0.5*magSqr(U))
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
    gradT_(fvc::grad(T_,"grad(TSlope)")),
    gradk_(fvc::grad(turbulenceModel_.k(),"grad(TKE)"))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux>
void Foam::godunovFlux<Flux>::update(Switch secondOrder)
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    const volVectorField& cellCenter = mesh_.C();
    const surfaceVectorField& faceCenter = mesh_.Cf();

    // need a copy here, because the return value is a <tmp>
//     const volScalarField Cv_ = thermophysicalModel_.Cv();
//     const volScalarField R_ =
//         ( thermophysicalModel_.Cp() - thermophysicalModel_.Cv() );
    const volScalarField kappa_ =
        thermophysicalModel_.Cp()/thermophysicalModel_.Cv();

    // only valid for ideal gas!
    volScalarField rho_("godunovRho",thermophysicalModel_.rho());
//     volScalarField rho_("godunovRho",p_/(R_*T_));
//     rho_.correctBoundaryConditions();
    volVectorField gradrho_ = fvc::grad(rho_,"grad(rhoSlope)");

    const volScalarField k_("TKE",turbulenceModel_.k());

    // Read field bounds
//     dictionary fieldBounds = mesh_.solutionDict().subDict("fieldBounds");

//     // Pressure bounds
//     dimensionedScalar pMin("pMin", p_.dimensions(), 0);
//     dimensionedScalar pMax("pMax", p_.dimensions(), 0);
//     fieldBounds.lookup(p_.name()) >> pMin.value() >> pMax.value();
// 
//     // Velocity bounds
//     dimensionedScalar smallU("smallU", dimVelocity, 1e-10);
//     dimensionedScalar UMax("UMax", U_.dimensions(), 0);
//     fieldBounds.lookup(U_.name()) >> UMax.value();
// 
//     // Temperature bounds
//     dimensionedScalar TMin("TMin", T_.dimensions(), 0);
//     dimensionedScalar TMax("TMax", T_.dimensions(), 0);
//     fieldBounds.lookup(T_.name()) >> TMin.value() >> TMax.value();

    // Switch between 1D and MultiDimensional Limiters
    bool multidimLimiterSwitch = false;

    // RTS 1D Limiters
    word limiterName("vanAlbadaSlope");

    // epsilon for the VK limiter
    word epsilon("5");
    // constant for Roe Entropy fix
    scalar Konstant = 0.05;

    // read riemann solver coeffs
    if(mesh_.solutionDict().found("Riemann"))
    {
        // Read Riemann coeffs dict
        dictionary riemann = mesh_.solutionDict().subDict("Riemann");
        if (riemann.found("multidimLimiter"))
        {
            multidimLimiterSwitch = Switch(riemann.lookup("multidimLimiter"));
        }
        if (riemann.found("limiterName"))
        {
            limiterName = word(riemann.lookup("limiterName"));
        }
        if (riemann.found("epsilon"))
        {
            epsilon = word(riemann.lookup("epsilon"));
        }
        Konstant = riemann.lookupOrDefault("RoeKonstant",Konstant);
    }

    // 2nd order correction
    IStringStream blendingFactor(epsilon);

    // MultiDimensional Limiters
//     BarthJespersenSlopeMultiLimiter scalarLimiter(blendingFactor);
//     BarthJespersenSlopeMultiLimiter vectorLimiter(blendingFactor);

    VenkatakrishnanSlopeMultiLimiter scalarLimiter(blendingFactor);
    VenkatakrishnanSlopeMultiLimiter vectorLimiter(blendingFactor);

    // 1D Limiters
//     BarthJespersenSlopeLimiter<NVDTVD> scalarLimiter(blendingFactor);
//     BarthJespersenSlopeLimiter<NVDVTVDV> vectorLimiter(blendingFactor);

//     MinmodSlopeLimiter<NVDTVD> scalarLimiter(blendingFactor);
//     MinmodSlopeLimiter<NVDVTVDV> vectorLimiter(blendingFactor);

//     vanAlbadaSlopeLimiter<NVDTVD> scalarLimiter(blendingFactor);
//     vanAlbadaSlopeLimiter<NVDVTVDV> vectorLimiter(blendingFactor);

//     vanLeerSlopeLimiter<NVDTVD> scalarLimiter(blendingFactor);
//     vanLeerSlopeLimiter<NVDVTVDV> vectorLimiter(blendingFactor);

    // Calculate and store the cell based limiter
    volScalarField pLimiter
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
    );

    volVectorField ULimiter
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
    );

//     volScalarField TLimiter
//     (
//         IOobject
//         (
//             "TLimiter",
//             mesh().time().timeName(),
//             mesh(),
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         mesh(),
//         dimensionedScalar("TLimiter", dimless, 1.0),
//         zeroGradientFvPatchScalarField::typeName
//     );

    volScalarField rhoLimiter
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
    );

    volScalarField kLimiter
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
    );

    if (secondOrder)
    {
        // Update the primitive gradients
        gradp_ = fvc::grad(p_,"grad(pSlope)");
        gradU_ = fvc::grad(U_,"grad(USlope)");
//         gradT_ = fvc::grad(T_,"grad(TSlope)");
        if (max(k_).value() > 0.0)
        {
            gradk_ = fvc::grad(k_);
        }

        if (multidimLimiterSwitch)
        {
            volScalarField pMinValue("pMin",p_);
            volScalarField pMaxValue("pMax",p_);
            volVectorField UMinValue("UMin",U_);
            volVectorField UMaxValue("UMax",U_);
//             volScalarField TMinValue("TMin",T_);
//             volScalarField TMaxValue("TMax",T_);

            volScalarField rhoMinValue("rhoMin",rho_);
            volScalarField rhoMaxValue("rhoMax",rho_);

            volScalarField kMinValue("kMin",k_);
            volScalarField kMaxValue("kMax",k_);

            // for BJ /VK find min/max of each value for each variable
            // internal face
            forAll(owner, faceI)
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                // min values
                pMinValue[own] =
                    min(pMinValue[own], p_[nei]);

                pMinValue[nei] =
                    min(pMinValue[nei], p_[own]);

                UMinValue[own] =
                    min(UMinValue[own], U_[nei]);

                UMinValue[nei] =
                    min(UMinValue[nei], U_[own]);

//                 TMinValue[own] =
//                     min(TMinValue[own], T_[nei]);
// 
//                 TMinValue[nei] =
//                     min(TMinValue[nei], T_[own]);

                rhoMinValue[own] =
                    min(rhoMinValue[own], rho_[nei]);

                rhoMinValue[nei] =
                    min(rhoMinValue[nei], rho_[own]);

                kMinValue[own] =
                    min(kMinValue[own], k_[nei]);

                kMinValue[nei] =
                    min(kMinValue[nei], k_[own]);

                // max values
                pMaxValue[own] =
                    max(pMaxValue[own], p_[nei]);

                pMaxValue[nei] =
                    max(pMaxValue[nei], p_[own]);

                UMaxValue[own] =
                    max(UMaxValue[own], U_[nei]);

                UMaxValue[nei] =
                    max(UMaxValue[nei], U_[own]);

//                 TMaxValue[own] =
//                     max(TMaxValue[own], T_[nei]);
// 
//                 TMaxValue[nei] =
//                     max(TMaxValue[nei], T_[own]);

                rhoMaxValue[own] =
                    max(rhoMaxValue[own], rho_[nei]);

                rhoMaxValue[nei] =
                    max(rhoMaxValue[nei], rho_[own]);

                kMaxValue[own] =
                    max(kMaxValue[own], k_[nei]);

                kMaxValue[nei] =
                    max(kMaxValue[nei], k_[own]);
            }

            // Update coupled boundary min/max values of primitive variables
            forAll(p_.boundaryField(), patchi)
            {
                const fvPatchScalarField& pp = p_.boundaryField()[patchi];
                const fvPatchVectorField& pU = U_.boundaryField()[patchi];
//                 const fvPatchScalarField& pT = T_.boundaryField()[patchi];
                const fvPatchScalarField& prho = rho_.boundaryField()[patchi];
                const fvPatchScalarField& pk = k_.boundaryField()[patchi];

                const unallocLabelList& faceCells = pp.patch().faceCells();

                if (pp.coupled())
                {
                    // primitive variables
                    const scalarField ppRight = pp.patchNeighbourField();
                    const vectorField pURight = pU.patchNeighbourField();
//                     const scalarField pTRight = pT.patchNeighbourField();
                    const scalarField prhoRight = prho.patchNeighbourField();
                    const scalarField pkRight = pk.patchNeighbourField();

                    forAll(pp, facei)
                    {
                        label own = faceCells[facei];

                        // min values at coupled boundary faces
                        pMinValue[own] =
                            min(pMinValue[own], ppRight[facei]);

                        UMinValue[own] =
                            min(UMinValue[own], pURight[facei]);

//                         TMinValue[own] =
//                             min(TMinValue[own], pTRight[facei]);

                        rhoMinValue[own] =
                            min(rhoMinValue[own], prhoRight[facei]);

                        kMinValue[own] =
                            min(kMinValue[own], pkRight[facei]);

                        // max values at coupled boundary faces
                        pMaxValue[own] =
                            max(pMaxValue[own], ppRight[facei]);

                        UMaxValue[own] =
                            max(UMaxValue[own], pURight[facei]);

//                         TMaxValue[own] =
//                             max(TMaxValue[own], pTRight[facei]);

                        rhoMaxValue[own] =
                            max(rhoMaxValue[own], prhoRight[facei]);

                        kMaxValue[own] =
                            max(kMaxValue[own], pkRight[facei]);
                    }
                }
            }

            // Not sure if this is really needed
            pMinValue.correctBoundaryConditions();
            UMinValue.correctBoundaryConditions();
//             TMinValue.correctBoundaryConditions();
            rhoMinValue.correctBoundaryConditions();
            kMinValue.correctBoundaryConditions();
            pMaxValue.correctBoundaryConditions();
            UMaxValue.correctBoundaryConditions();
//             TMaxValue.correctBoundaryConditions();
            rhoMaxValue.correctBoundaryConditions();
            kMaxValue.correctBoundaryConditions();

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

            // compute for each cell a limiter
            // Loop over all faces, which very inefficent, as the limiter is
            // computed several times!
            forAll(owner, faceI)
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                vector deltaRLeft  = faceCenter[faceI] - cellCenter[own];
                vector deltaRRight = faceCenter[faceI] - cellCenter[nei];

                // flux dummy
                scalar upwindFlux = 1.0;

                scalar pOwnerLimiter = scalarLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,
                    pMaxValue[own]-p_[own],
                    pMinValue[own]-p_[own],
                    gradp_[own],
                    gradp_[own],
                    deltaRLeft
                );
                scalar pNeighbourLimiter = scalarLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    pMaxValue[nei]-p_[nei],
                    pMinValue[nei]-p_[nei],
                    gradp_[nei],
                    gradp_[nei],
                    deltaRRight
                );

                vector UOwnerLimiterV = vectorLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,
                    UMaxValue[own]-U_[own],
                    UMinValue[own]-U_[own],
                    gradU_[own],
                    gradU_[own],
                    deltaRLeft
                );
                vector UNeighbourLimiterV = vectorLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    UMaxValue[nei]-U_[nei],
                    UMinValue[nei]-U_[nei],
                    gradU_[nei],
                    gradU_[nei],
                    deltaRRight
                );

//                 scalar TOwnerLimiter = scalarLimiter.limiter
//                 (
//                     cellVolume[own],
//                     upwindFlux,
//                     TMaxValue[own]-T_[own],
//                     TMinValue[own]-T_[own],
//                     gradT_[own],
//                     gradT_[own],
//                     deltaRLeft
//                 );
//                 scalar TNeighbourLimiter = scalarLimiter.limiter
//                 (
//                     cellVolume[nei],
//                     upwindFlux,
//                     TMaxValue[nei]-T_[nei],
//                     TMinValue[nei]-T_[nei],
//                     gradT_[nei],
//                     gradT_[nei],
//                     deltaRRight
//                 );

                scalar rhoOwnerLimiter = scalarLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,
                    rhoMaxValue[own]-rho_[own],
                    rhoMinValue[own]-rho_[own],
                    gradrho_[own],
                    gradrho_[own],
                    deltaRLeft
                );
                scalar rhoNeighbourLimiter = scalarLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    rhoMaxValue[nei]-rho_[nei],
                    rhoMinValue[nei]-rho_[nei],
                    gradrho_[nei],
                    gradrho_[nei],
                    deltaRRight
                );

                scalar kOwnerLimiter = scalarLimiter.limiter
                (
                    cellVolume[own],
                    upwindFlux,
                    kMaxValue[own]-k_[own],
                    kMinValue[own]-k_[own],
                    gradk_[own],
                    gradk_[own],
                    deltaRLeft
                );
                scalar kNeighbourLimiter = scalarLimiter.limiter
                (
                    cellVolume[nei],
                    upwindFlux,
                    kMaxValue[nei]-k_[nei],
                    kMinValue[nei]-k_[nei],
                    gradk_[nei],
                    gradk_[nei],
                    deltaRRight
                );

                // find minimal limiter value in each cell
                pLimiter[own] =
                    min(pLimiter[own], pOwnerLimiter);

                pLimiter[nei] =
                    min(pLimiter[nei], pNeighbourLimiter);

                ULimiter[own] =
                    min(ULimiter[own], UOwnerLimiterV);

                ULimiter[nei] =
                    min(ULimiter[nei], UNeighbourLimiterV);

//                 TLimiter[own] =
//                     min(TLimiter[own], TOwnerLimiter);
// 
//                 TLimiter[nei] =
//                     min(TLimiter[nei], TNeighbourLimiter);

                rhoLimiter[own] =
                    min(rhoLimiter[own], rhoOwnerLimiter);

                rhoLimiter[nei] =
                    min(rhoLimiter[nei], rhoNeighbourLimiter);

                kLimiter[own] =
                    min(kLimiter[own], kOwnerLimiter);

                kLimiter[nei] =
                    min(kLimiter[nei], kNeighbourLimiter);
            }

            // Update coupled boundary limiters
            forAll(p_.boundaryField(), patchi)
            {
                const fvPatchScalarField& pp = p_.boundaryField()[patchi];

//                 const fvPatchVectorField& pCellCenter = cellCenter.boundaryField()[patchi];
                const unallocLabelList& faceCells = pp.patch().faceCells();

                if (pp.coupled())
                {
                    // cell and face centers
                    const vectorField delta =  pp.patch().delta();
//                     const vectorField faceCenter =  pp.patch().Cf();
//                     const vectorField pCellCenterLeft  =  pCellCenter.patchInternalField();

                    forAll(pp, facei)
                    {
                        label own = faceCells[facei];

//                     vector deltaRLeft  = faceCenter[facei] - pCellCenterLeft[facei];
                        vector deltaRLeft  = delta[facei];

                        // flux dummy
                        scalar upwindFlux = 1.0;

                        scalar pOwnerLimiter = scalarLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            pMaxValue[own]-p_[own],
                            pMinValue[own]-p_[own],
                            gradp_[own],
                            gradp_[own],
                            deltaRLeft
                        );
                        vector UOwnerLimiterV = vectorLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            UMaxValue[own]-U_[own],
                            UMinValue[own]-U_[own],
                            gradU_[own],
                            gradU_[own],
                            deltaRLeft
                        );
//                         scalar TOwnerLimiter = scalarLimiter.limiter
//                         (
//                             cellVolume[own],
//                             upwindFlux,
//                             TMaxValue[own]-T_[own],
//                             TMinValue[own]-T_[own],
//                             gradT_[own],
//                             gradT_[own],
//                             deltaRLeft
//                         );
                        scalar rhoOwnerLimiter = scalarLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            rhoMaxValue[own]-rho_[own],
                            rhoMinValue[own]-rho_[own],
                            gradrho_[own],
                            gradrho_[own],
                            deltaRLeft
                        );
                        scalar kOwnerLimiter = scalarLimiter.limiter
                        (
                            cellVolume[own],
                            upwindFlux,
                            kMaxValue[own]-k_[own],
                            kMinValue[own]-k_[own],
                            gradk_[own],
                            gradk_[own],
                            deltaRLeft
                        );

                        pLimiter[own] =
                            min(pLimiter[own], pOwnerLimiter);
                        ULimiter[own] =
                            min(ULimiter[own], UOwnerLimiterV);
//                         TLimiter[own] =
//                             min(TLimiter[own], TOwnerLimiter);
                        rhoLimiter[own] =
                            min(rhoLimiter[own], rhoOwnerLimiter);
                        kLimiter[own] =
                            min(kLimiter[own], kOwnerLimiter);
                    }
                }
//                 else
//                 {
//                     forAll(pp, facei)
//                     {
//                         label own = faceCells[facei];
//                         // Calculate minimal limiters
//                         pLimiter[own] =
//                             min(pLimiter[own], 1);
//                         ULimiter[faceCells[facei]] =
//                             min(ULimiter[own], vector::one);
//                         TLimiter[faceCells[facei]] =
//                             min(TLimiter[own], 1);
//                     }
//                 }
            }
        }
        else
        {
            IStringStream pSchemeData(limiterName);
            IStringStream USchemeData(limiterName+"V");
//             IStringStream TSchemeData(limiterName);
            IStringStream rhoSchemeData(limiterName);
            IStringStream kSchemeData(limiterName);

            // is working
            surfaceScalarField pLimiterSurface =
                limitedSurfaceInterpolationScheme<scalar>::New
                (
                    mesh_,
                    rhoFlux_,
                    pSchemeData
                )().limiter(p_);

            surfaceScalarField ULimiterSurface =
                limitedSurfaceInterpolationScheme<vector>::New
                (
                    mesh_,
                    rhoFlux_,
                    USchemeData
                )().limiter(U_);

//             surfaceScalarField TLimiterSurface =
//                 limitedSurfaceInterpolationScheme<scalar>::New
//                 (
//                     mesh_,
//                     rhoFlux_,
//                     TSchemeData
//                 )().limiter(T_);

            surfaceScalarField rhoLimiterSurface =
                limitedSurfaceInterpolationScheme<scalar>::New
                (
                    mesh_,
                    rhoFlux_,
                    rhoSchemeData
                )().limiter(rho_);

            surfaceScalarField kLimiterSurface =
                limitedSurfaceInterpolationScheme<scalar>::New
                (
                    mesh_,
                    rhoFlux_,
                    kSchemeData
                )().limiter(k_);

            // 1D implementation
            forAll(owner, faceI)
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                // find minimal limiter value in each cell
                pLimiter[own] =
                    min(pLimiter[own], pLimiterSurface[faceI]);

                pLimiter[nei] =
                    min(pLimiter[nei], pLimiterSurface[faceI]);

                ULimiter[own] =
                    min(ULimiter[own], ULimiterSurface[faceI]*vector::one);

                ULimiter[nei] =
                    min(ULimiter[nei], ULimiterSurface[faceI]*vector::one);

//                 TLimiter[own] =
//                     min(TLimiter[own], TLimiterSurface[faceI]);
// 
//                 TLimiter[nei] =
//                     min(TLimiter[nei], TLimiterSurface[faceI]);

                rhoLimiter[own] =
                    min(rhoLimiter[own], rhoLimiterSurface[faceI]);

                rhoLimiter[nei] =
                    min(rhoLimiter[nei], rhoLimiterSurface[faceI]);

                kLimiter[own] =
                    min(kLimiter[own], kLimiterSurface[faceI]);

                kLimiter[nei] =
                    min(kLimiter[nei], kLimiterSurface[faceI]);
            }

            // Update coupled boundary limiters
            forAll(p_.boundaryField(), patchi)
            {
                const fvPatchScalarField& pp = p_.boundaryField()[patchi];
                const unallocLabelList& faceCells = pp.patch().faceCells();

                if (pp.coupled())
                {
                    forAll(pp, facei)
                    {
                        label own = faceCells[facei];

                        pLimiter[own] =
                            min(pLimiter[own], pLimiterSurface[facei]);
                        ULimiter[own] =
                            min(ULimiter[own], ULimiterSurface[facei]*vector::one);
//                         TLimiter[own] =
//                             min(TLimiter[own], TLimiterSurface[facei]);
                        rhoLimiter[own] =
                            min(rhoLimiter[own], rhoLimiterSurface[facei]);
                        kLimiter[own] =
                            min(kLimiter[own], kLimiterSurface[facei]);
                    }
                }
            }
        }
    }

    // Not sure if this is really needed
    pLimiter.correctBoundaryConditions();
    ULimiter.correctBoundaryConditions();
//     TLimiter.correctBoundaryConditions();
    rhoLimiter.correctBoundaryConditions();
    kLimiter.correctBoundaryConditions();

//     boundMinMax(pLimiter,dimensionedScalar(0.0),dimensionedScalar(1.0));
//     boundMinMax(ULimiter[0],dimensionedScalar(0.0),dimensionedScalar(1.0));
//     boundMinMax(TLimiter,dimensionedScalar(0.0),dimensionedScalar(1.0));

    Info << "max(pLimiter) "  << max(pLimiter.internalField())
         << " min(pLimiter) " << min(pLimiter.internalField())  << endl;
    Info << "max(ULimiter) "  << max(ULimiter.internalField())
         << " min(ULimiter) " << min(ULimiter.internalField())  << endl;
//     Info << "max(TLimiter) "  << max(TLimiter.internalField())
//          << " min(TLimiter) " << min(TLimiter.internalField())  << endl;
    Info << "max(rhoLimiter) "  << max(rhoLimiter.internalField())
         << " min(rhoLimiter) " << min(rhoLimiter.internalField())  << endl;
    Info << "max(kLimiter) "  << max(kLimiter.internalField())
         << " min(kLimiter) " << min(kLimiter.internalField())  << endl;

    // Calculate fluxes at internal faces
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector deltaRLeft  = faceCenter[faceI] - cellCenter[own];
        vector deltaRRight = faceCenter[faceI] - cellCenter[nei];

        // calculate fluxes with reconstructed primitive variables at faces
        // TODO: thermophysical variables are not reconstructed at faces!!!
        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
//             p_[own],      // left p
//             p_[nei],      // right p
//             U_[own],      // left U
//             U_[nei],      // right U
//             T_[own],      // left T
//             T_[nei],      // right T
//
//             max(p_[own] + secondOrder*pLimiter[own]*(gradp_[own] & deltaRLeft),pMin.value()), // reconstructed left p
//             max(p_[nei] + secondOrder*pLimiter[nei]*(gradp_[nei] & deltaRRight),pMin.value()),// reconstructed right p
//
            p_[own] + secondOrder*pLimiter[own]*(gradp_[own] & deltaRLeft), // reconstructed left p
            p_[nei] + secondOrder*pLimiter[nei]*(gradp_[nei] & deltaRRight),// reconstructed right p
//
//             U_[own] + pLimiter[own] * (gradU_[own] & deltaRLeft), // reconstructed left U
//             U_[nei] + pLimiter[nei] * (gradU_[nei] & deltaRRight),// reconstructed right U
//             T_[own] + pLimiter[own] * (gradT_[own] & deltaRLeft), // reconstructed left T
//             T_[nei] + pLimiter[nei] * (gradT_[nei] & deltaRRight),// reconstructed right T
//
            U_[own] + secondOrder*cmptMultiply(ULimiter[own],(gradU_[own] & deltaRLeft)), // reconstructed left U
            U_[nei] + secondOrder*cmptMultiply(ULimiter[nei],(gradU_[nei] & deltaRRight)),// reconstructed right U
//
//          using minimum component scalar limiter
//             U_[own] + secondOrder*cmptMin(ULimiter[own])*(gradU_[own] & deltaRLeft), // reconstructed left U
//             U_[nei] + secondOrder*cmptMin(ULimiter[nei])*(gradU_[nei] & deltaRRight),// reconstructed right U
//
//             U_[own] + secondOrder*ULimiter[own]*(gradU_[own] & deltaRLeft), // reconstructed left U
//             U_[nei] + secondOrder*ULimiter[nei]*(gradU_[nei] & deltaRRight),// reconstructed right U
//
//             max(T_[own] + secondOrder*TLimiter[own]*(gradT_[own] & deltaRLeft),TMin.value()), // reconstructed left T
//             max(T_[nei] + secondOrder*TLimiter[nei]*(gradT_[nei] & deltaRRight),TMin.value()),// reconstructed right T
//
//             T_[own] + secondOrder*TLimiter[own]*(gradT_[own] & deltaRLeft), // reconstructed left T
//             T_[nei] + secondOrder*TLimiter[nei]*(gradT_[nei] & deltaRRight),// reconstructed right T
//
            rho_[own] + secondOrder*rhoLimiter[own]*(gradrho_[own] & deltaRLeft), // reconstructed left T
            rho_[nei] + secondOrder*rhoLimiter[nei]*(gradrho_[nei] & deltaRRight),// reconstructed right T
//
            k_[own] + secondOrder*kLimiter[own]*(gradk_[own] & deltaRLeft),
            k_[nei] + secondOrder*kLimiter[nei]*(gradk_[nei] & deltaRRight),
// //
//             R_[own],        // left R
//             R_[nei],        // right R
//             Cv_[own],       // left Cv
//             Cv_[nei],       // right Cv
// //
            kappa_[own],       // left kappa
            kappa_[nei],       // right kappa
// //
            Sf[faceI],      // face vector
            magSf[faceI],   // face area
            dotX_[faceI],    // face velocity
            Konstant
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
//         const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        const fvPatchScalarField& prho = rho_.boundaryField()[patchi];
        const fvPatchScalarField& pk = k_.boundaryField()[patchi];

        const fvPatchVectorField& pGradp = gradp_.boundaryField()[patchi];
        const fvPatchTensorField& pGradU = gradU_.boundaryField()[patchi];
//         const fvPatchVectorField& pGradT = gradT_.boundaryField()[patchi];
        const fvPatchVectorField& pGradrho = gradrho_.boundaryField()[patchi];
        const fvPatchVectorField& pGradk = gradk_.boundaryField()[patchi];

//         const fvPatchScalarField& pCv = Cv_.boundaryField()[patchi];
//         const fvPatchScalarField& pR = R_.boundaryField()[patchi];
        const fvPatchScalarField& pkappa = kappa_.boundaryField()[patchi];

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const fvsPatchVectorField& pDotX = dotX_.boundaryField()[patchi];

        const fvPatchVectorField& pCellCenter = cellCenter.boundaryField()[patchi];
//         const fvsPatchVectorField& pFaceCenter = faceCenter.boundaryField()[patchi];

        const fvPatchScalarField& ppLimiter = pLimiter.boundaryField()[patchi];
        const fvPatchVectorField& pULimiter = ULimiter.boundaryField()[patchi];
//         const fvPatchScalarField& pTLimiter = TLimiter.boundaryField()[patchi];
        const fvPatchScalarField& prhoLimiter = rhoLimiter.boundaryField()[patchi];
        const fvPatchScalarField& pkLimiter = kLimiter.boundaryField()[patchi];

        if (pp.coupled())
        {
            // primitive variables
            const scalarField ppLeft  = pp.patchInternalField();
            const scalarField ppRight = pp.patchNeighbourField();

            const vectorField pULeft  = pU.patchInternalField();
            const vectorField pURight = pU.patchNeighbourField();

//             const scalarField pTLeft  = pT.patchInternalField();
//             const scalarField pTRight = pT.patchNeighbourField();

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

//             const vectorField pGradTLeft  = pGradT.patchInternalField();
//             const vectorField pGradTRight = pGradT.patchNeighbourField();

            const vectorField pGradrhoLeft  = pGradrho.patchInternalField();
            const vectorField pGradrhoRight = pGradrho.patchNeighbourField();

            const vectorField pGradkLeft  = pGradk.patchInternalField();
            const vectorField pGradkRight = pGradk.patchNeighbourField();

            // cell limiters
            scalarField ppLimiterLeft  = ppLimiter.patchInternalField();
            scalarField ppLimiterRight = ppLimiter.patchNeighbourField();

            vectorField pULimiterLeft  = pULimiter.patchInternalField();
            vectorField pULimiterRight = pULimiter.patchNeighbourField();

//             scalarField pTLimiterLeft  = pTLimiter.patchInternalField();
//             scalarField pTLimiterRight = pTLimiter.patchNeighbourField();

            scalarField prhoLimiterLeft  = prhoLimiter.patchInternalField();
            scalarField prhoLimiterRight = prhoLimiter.patchNeighbourField();

            scalarField pkLimiterLeft  = pkLimiter.patchInternalField();
            scalarField pkLimiterRight = pkLimiter.patchNeighbourField();

            // cell and face centers
            const vectorField faceCenter =  pp.patch().Cf();
//             const vectorField deltaRLeftField = pp.patch().delta();
            const vectorField pCellCenterLeft  =  pCellCenter.patchInternalField();
            const vectorField pCellCenterRight =  pCellCenter.patchNeighbourField();

            forAll(pp, facei)
            {
                vector deltaRLeft  = faceCenter[facei] - pCellCenterLeft[facei];
                vector deltaRRight = faceCenter[facei] - pCellCenterRight[facei];

                // bound the limiters between 0 and 1 in order to prevent problems due to interpolations

                ppLimiterLeft[facei]   = max(min(ppLimiterLeft[facei],1.0),0.0);
                ppLimiterRight[facei]  = max(min(ppLimiterRight[facei],1.0),0.0);

                pULimiterLeft[facei]  = max(min(pULimiterLeft[facei],vector::one),vector::zero);
                pULimiterRight[facei] = max(min(pULimiterRight[facei],vector::one),vector::zero);

//                 pTLimiterLeft[facei]   = max(min(pTLimiterLeft[facei],1.0),0.0);
//                 pTLimiterRight[facei]  = max(min(pTLimiterRight[facei],1.0),0.0);

                prhoLimiterLeft[facei]   = max(min(prhoLimiterLeft[facei],1.0),0.0);
                prhoLimiterRight[facei]  = max(min(prhoLimiterRight[facei],1.0),0.0);

                pkLimiterLeft[facei]   = max(min(pkLimiterLeft[facei],1.0),0.0);
                pkLimiterRight[facei]  = max(min(pkLimiterRight[facei],1.0),0.0);

                // Calculate fluxes at coupled boundary faces
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
//
//                     ppLeft[facei],    // face p
//                     ppRight[facei],   // face p
//                     pULeft[facei],    // face U
//                     pURight[facei],   // face U
//                     pTLeft[facei],    // face T
//                     pTRight[facei],   // face T
//                     pkLeft[facei],    // face k
//                     pkRight[facei],   // face k
//
// //                     max(ppLeft[facei]  + secondOrder*ppLimiterLeft[facei]*(pGradpLeft[facei] & deltaRLeft),pMin.value()),                  // face p
// //                     max(ppRight[facei] + secondOrder*ppLimiterRight[facei]*(pGradpRight[facei] & deltaRRight),pMin.value()),               // face p
//
                    ppLeft[facei]  + secondOrder*ppLimiterLeft[facei]*(pGradpLeft[facei] & deltaRLeft),                  // face p
                    ppRight[facei] + secondOrder*ppLimiterRight[facei]*(pGradpRight[facei] & deltaRRight),               // face p
// //
                    pULeft[facei]  + secondOrder*cmptMultiply(pULimiterLeft[facei],(pGradULeft[facei] & deltaRLeft)),    // face U
                    pURight[facei] + secondOrder*cmptMultiply(pULimiterRight[facei],(pGradURight[facei] & deltaRRight)), // face U
// //
//                     pULeft[facei]  + secondOrder*cmptMin(pULimiterLeft[facei])*(pGradULeft[facei] & deltaRLeft),    // face U
//                     pURight[facei] + secondOrder*cmptMin(pULimiterRight[facei])*(pGradURight[facei] & deltaRRight), // face U
// //
//                     max(pTLeft[facei]  + secondOrder*pTLimiterLeft[facei]*(pGradTLeft[facei] & deltaRLeft),TMin.value()),                 // face T
//                     max(pTRight[facei] + secondOrder*pTLimiterRight[facei]*(pGradTRight[facei] & deltaRRight),TMin.value()),               // face T
//
//                     pTLeft[facei]  + secondOrder*pTLimiterLeft[facei] *(pGradTLeft[facei]  & deltaRLeft),                 // face T
//                     pTRight[facei] + secondOrder*pTLimiterRight[facei]*(pGradTRight[facei] & deltaRRight),               // face T
//
                    prhoLeft[facei]  + secondOrder*prhoLimiterLeft[facei] *(pGradrhoLeft[facei]  & deltaRLeft),                 // face T
                    prhoRight[facei] + secondOrder*prhoLimiterRight[facei]*(pGradrhoRight[facei] & deltaRRight),               // face T
//
                    pkLeft[facei]  + secondOrder*pkLimiterLeft[facei]*(pGradkLeft[facei] & deltaRLeft),                 // face T
                    pkRight[facei] + secondOrder*pkLimiterRight[facei]*(pGradkRight[facei] & deltaRRight),               // face T
// //
//                     pR[facei],        // face R
//                     pR[facei],        // face R
//                     pCv[facei],       // face Cv
//                     pCv[facei],       // face Cv
// //
                    pkappaLeft[facei],       // face kappa
                    pkappaRight[facei],       // face kappa
// //
                    pSf[facei],       // face vector
                    pMagSf[facei],    // face area
                    pDotX[facei],      // face velocity
                    Konstant
                );
            }
        }
        else
        {
            forAll(pp, facei)
            {
                // Calculate fluxes at boundary faces
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pp[facei],        // face p
                    pp[facei],        // face p
                    pU[facei],        // face U
                    pU[facei],        // face U
// //
//                     pT[facei],        // face T
//                     pT[facei],        // face T
// //
                    prho[facei],        // face T
                    prho[facei],        // face T
// //
                    pk[facei],        // face k
                    pk[facei],        // face k
// //
//                     pR[facei],        // face R
//                     pR[facei],        // face R
//                     pCv[facei],       // face Cv
//                     pCv[facei],       // face Cv
// //
                    pkappa[facei],       // face kappa
                    pkappa[facei],       // face kappa
// //
                    pSf[facei],       // face vector
                    pMagSf[facei],    // face area
                    pDotX[facei],      // face velocity
                    Konstant
                );
            }
        }
    }
}

// ************************************************************************* //
