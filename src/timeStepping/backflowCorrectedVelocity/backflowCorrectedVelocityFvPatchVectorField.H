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
    Foam::backflowCorrectedVelocityFvPatchVectorField

Description
    Velocity outlet boundary condition for patches where the pressure is
    specified. The outflow velocity is obtained by "zeroGradient".
    Reverse flow is cutted with the preFactor variable. If preFactor is larger
    one, positive flow is applied.
    Should only be used in conjunction with outletInlet or outletInletBackflow
    boundary condition for the pressure.

SourceFiles
    backflowCorrectedVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef backflowCorrectedVelocityFvPatchVectorField_H
#define backflowCorrectedVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fluxCorrectedVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class backflowCorrectedVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
    // Private data

        //- Scalar preFactor
        scalar preFactor_;

public:

    //- Runtime type information
    TypeName("backflowCorrectedVelocity");


    // Constructors

        //- Construct from patch and internal field
        backflowCorrectedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        backflowCorrectedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given backflowCorrectedVelocityFvPatchVectorField
        //  onto a new patch
        backflowCorrectedVelocityFvPatchVectorField
        (
            const backflowCorrectedVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new backflowCorrectedVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        backflowCorrectedVelocityFvPatchVectorField
        (
            const backflowCorrectedVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new backflowCorrectedVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
