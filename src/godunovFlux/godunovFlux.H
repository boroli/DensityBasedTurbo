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
    Generic Godunov flux class. Does an extrapolation of cell centered
    variables to face centeres, using a Taylor-Series Expansion. Therefore a
    slope limiter is utilized.

    TODO: Make the slope limiter run-time selectable.

    References:
    Weiss, Jonathan M. and Wayne A. Smith, "Preconditioning Applied to Variable
    and Constant Density Flows." AIAA Journal 33/11 (November 1995): 2050–2057.

Author
    Oliver Borm  All rights reserved.

SourceFiles
    godunovFlux.H
    godunovFlux.C

\*---------------------------------------------------------------------------*/

#ifndef godunovFlux_H
#define godunovFlux_H

#include "fvCFD.H"
#include "basicThermo.H"
#include "turbulenceModel.H"

#include "SlopeLimiter/BarthJespersenSlopeMulti/BarthJespersenSlopeMulti.H"
#include "SlopeLimiter/VenkatakrishnanSlopeMulti/VenkatakrishnanSlopeMulti.H"

// #include "hllcFlux/hllcFlux.H"
#include "hllcALEFlux/hllcALEFlux.H"
// #include "hllcALEDGFlux/hllcALEDGFlux.H"
// #include "roeFlux/roeFlux.H"
#include "roeALEFlux/roeALEFlux.H"
#include "rusanovALEFlux/rusanovALEFlux.H"
#include "AUSMplusALEFlux/AUSMplusALEFlux.H"
#include "AUSMplusUpALEFlux/AUSMplusUpALEFlux.H"

#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class godunovFlux Declaration
\*---------------------------------------------------------------------------*/

template<class Flux>
class godunovFlux
:
    public Flux
{
    // Private data

        //- Reference to a flux class
//         const Flux& flux_;

        //- Reference to mesh
        const fvMesh& mesh_;

        // Reference to primitive fields

            //- static Pressure
            const volScalarField& p_;

            //- Velocity
            const  volVectorField& U_;

            //- static density
            const  volScalarField& rho_;

            //- Reference to the thermophysicalModel
            const basicThermo& thermophysicalModel_;

            //- Reference to the turbulenceModel
            const compressible::turbulenceModel& turbulenceModel_;

            //- Get face-to-cell addressing: face area points from owner to neighbour
//             const unallocLabelList& owner_;

            //- Get face-to-cell addressing: face area points from owner to neighbour
//             const unallocLabelList& neighbour_;

            //- Reference to the cell centres
//             const volVectorField& cellCenter_;

            //- Reference to the face centres
//             const surfaceVectorField& faceCenter_;

            //- local copy of turbulent kinetic energy
//             volScalarField k_;

            //- local copy of isentropic expansion factor
            volScalarField kappa_;

            //- store cell volume
            volScalarField cellVolume_;

        // Fluxes

            //- Density flux
            surfaceScalarField rhoFlux_;

            //- Velocity flux
            surfaceVectorField rhoUFlux_;

            //- Energy flux
            surfaceScalarField rhoEFlux_;

            //- mesh fluxes
            surfaceVectorField dotX_;

        // Gradients

            //- static Pressure
            volVectorField gradp_;

            //- Velocity
            volTensorField gradU_;

            //- static density
            volVectorField gradrho_;

            //- turbulent kinetic energy
            volVectorField gradk_;

        // Limiters

            //- field for p Limiter
            volScalarField pLimiter_;

            //- field for U Limiter
            volVectorField ULimiter_;

            //- field for rho Limiter
            volScalarField rhoLimiter_;

            //- field for k Limiter
            volScalarField kLimiter_;

        // Parameters

            //- run-time-selectable parameter for the VK limiter
            scalar epsilon_;

            //- run-time-selectable parameter for the Riemann-Solvers:
            // Roe-Riemann-Solver: Konstant controls Harten's Entropy-Fix parameter
            // AUSM+-up-Riemann-Solver: Konstant sets Ma_infinity (MaInf)
            scalar Konstant_;

            //- 1D Limiter name
            word limiterName_;

            //- switch between 1D and multiDimensional limiter
            bool multidimLimiterSwitch_;

    // Private Member Functions

        //- compute Limiter field
        template<class phiType, class gradPhiType, class Limiter>
        void updateLimiter
        (
            const GeometricField<phiType, fvPatchField, volMesh>& phi,
            const GeometricField<gradPhiType, fvPatchField, volMesh>& gradPhi,
            GeometricField<phiType, fvPatchField, volMesh>& phiLimiter,
            word oneDLimiterName
        );

        //- Disallow default bitwise copy construct
        godunovFlux(const godunovFlux&);

        //- Disallow default bitwise assignment
        void operator=(const godunovFlux&);

public:

    // Constructors

        //- Construct from components
        godunovFlux
        (
            const volScalarField& p,
            const volVectorField& U,
            const volScalarField& rho,
            const basicThermo& thermo,
            const compressible::turbulenceModel& turbulenceModel
        );


    // Destructor - default

    // Member Functions

        //- Return reference to the mesh
        const fvMesh& mesh() const
        {
             return mesh_;
        }

        // Return fluxes

            //- Return density flux
            const surfaceScalarField& rhoFlux() const
            {
                return rhoFlux_;
            }

            //- Return velocity flux
            const surfaceVectorField& rhoUFlux() const
            {
                return rhoUFlux_;
            }

            //- Return energy flux
            const surfaceScalarField& rhoEFlux() const
            {
                return rhoEFlux_;
            }

            //- Return mesh flux
            const surfaceVectorField& dotX() const
            {
                return dotX_;
            }

            //- Return mesh flux
            surfaceVectorField& dotX()
            {
                return dotX_;
            }

        // Return Gradients

            //- Return pressure gradient
            const volVectorField& gradp() const
            {
                return gradp_;
            }

            //- Return pressure gradient
            const volTensorField& gradU() const
            {
                return gradU_;
            }

            //- Return TKE gradient
            const volVectorField& gradk() const
            {
                return gradk_;
            }

        //- Update fluxes based on current state
        void update(Switch secondOrder);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "godunovFlux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
