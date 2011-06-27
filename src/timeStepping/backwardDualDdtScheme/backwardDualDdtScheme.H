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
    Foam::fv::backwardDualDdtScheme

Description
    This is a dual-time stepping ddt scheme. For the physical time integration
    a second-order backward-differencing scheme is used. Whereas for pseudo
    time a forward Euler scheme with local time-stepping is used. This is
    cooked up with the current pseudo-time step and two previous physical
    time-step values. The following integration procedure applies for any
    GeoField Q, with the pseudo-time (subCycle) index i and the physical time
    index n.

    \f$  Q^{i+1} \left( \frac{1}{\beta_i \Delta \tau}
         + \frac{3}{2} \frac{1}{\Delta t} \right)
         - \left[ \frac{1}{\beta_i \Delta \tau} Q^{i}
         + \frac{4}{2} \frac{1}{\Delta t} Q^n
         - \frac{1}{2} \frac{1}{\Delta t} Q^{n-1} \right]
         = - R^i
    \f$

    These two previous physical time-step has to be allocated
    and prepared outside this ddt scheme! Also the local-time step field as
    well as the physical time steps and the multi-stage coefficient has to be
    provided and allocated outside of this ddt scheme! You have been warned!

    If the multi-stage coefficient is set to 1, this scheme reduces to a one
    step forward Euler scheme.

    TODO: Up to now only the three fvm::ddt functions are rewritten, all other
    methods are just copied from the standard Foam::fv::backwardDdtScheme

    TODO: Check the old cell volumes from the previous physical time steps for
    a moving mesh computation. And also the implementation of the mesh.phi()!

    References:
    Weiss, Jonathan M. and Wayne A. Smith, "Preconditioning Applied to Variable
    and Constant Density Flows." AIAA Journal 33/11 (November 1995): 2050–2057.

    Jameson, A. "Time Dependent Calculations Using Multigrid, with Applications
    to Unsteady Flows Past Airfoils and Wings" AIAA 91-1596

    Arnone, Liou and Povinelli, "Multigrid Time-Accurate Integration of
    Navier-Stokes Equations" AIAA 1993-3361, NASA-TM-106373

Author
    Oliver Borm  All rights reserved.

SourceFiles
    backwardDualDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef backwardDualDdtScheme_H
#define backwardDualDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class backwardDualDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class backwardDualDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Name of the physical time-step (field)
        word deltaTName_;

        //- Name of the pseudo time-step field
        word deltaTauName_;

        //- Name of old field
        word oldName_;

        //- current physical time-step
        scalar physDeltaT_;

        //- previous physical time-step
        scalar physDeltaT0_;

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the current physical time-step
        template<class GeoField>
        scalar deltaT_(const GeoField&) const;

        //- Return the previous physical time-step
        template<class GeoField>
        scalar deltaT0_(const GeoField&, const GeoField&) const;

//         scalar deltaT0_(IOList<label>) const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField&) const;

        //- Disallow default bitwise copy construct
        backwardDualDdtScheme(const backwardDualDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const backwardDualDdtScheme&);


public:

    //- Runtime type information
    TypeName("backwardDual");


    // Constructors

        //- Construct from mesh
        backwardDualDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh),
            deltaTName_("physDeltaT"),
            deltaTauName_("CoDeltaT"),
            oldName_("Old"),
            physDeltaT_(0.1),
            physDeltaT0_(0.1)
        {}

        //- Construct from mesh and Istream
        backwardDualDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is),
            deltaTName_(is),
            deltaTauName_(is),
            oldName_("Old"),
            physDeltaT_(0.1),
            physDeltaT0_(0.1)
        {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const dimensioned<Type>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            const dimensionedScalar&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            const volScalarField&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rA,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rA,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );


        tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


template<>
tmp<surfaceScalarField> backwardDualDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& U,
    const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> backwardDualDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "backwardDualDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
