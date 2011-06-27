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

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "temperatureDirectedInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

temperatureDirectedInletVelocityFvPatchVectorField::
temperatureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    TName_("T"),
    T0_(p.size()),
    inletDir_(p.size()),
    cylindricalCS_(0),
    omega_(vector::zero)
{}


temperatureDirectedInletVelocityFvPatchVectorField::
temperatureDirectedInletVelocityFvPatchVectorField
(
    const temperatureDirectedInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    TName_(ptf.TName_),
    T0_(ptf.T0_, mapper),
    inletDir_(ptf.inletDir_, mapper),
    cylindricalCS_(ptf.cylindricalCS_),
    omega_(ptf.omega_)
{}


temperatureDirectedInletVelocityFvPatchVectorField::
temperatureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    T0_("T0", dict, p.size()),
    inletDir_("inletDirection", dict, p.size()),
    cylindricalCS_(dict.lookup("cylindricalCS")),
    omega_(dict.lookup("omega"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


temperatureDirectedInletVelocityFvPatchVectorField::
temperatureDirectedInletVelocityFvPatchVectorField
(
    const temperatureDirectedInletVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    TName_(pivpvf.TName_),
    T0_(pivpvf.T0_),
    inletDir_(pivpvf.inletDir_),
    cylindricalCS_(pivpvf.cylindricalCS_),
    omega_(pivpvf.omega_)
{}


temperatureDirectedInletVelocityFvPatchVectorField::
temperatureDirectedInletVelocityFvPatchVectorField
(
    const temperatureDirectedInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    TName_(pivpvf.TName_),
    T0_(pivpvf.T0_),
    inletDir_(pivpvf.inletDir_),
    cylindricalCS_(pivpvf.cylindricalCS_),
    omega_(pivpvf.omega_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void temperatureDirectedInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    T0_.autoMap(m);
    inletDir_.autoMap(m);
}


void temperatureDirectedInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const temperatureDirectedInletVelocityFvPatchVectorField& tiptf =
        refCast<const temperatureDirectedInletVelocityFvPatchVectorField>(ptf);

    T0_.rmap(tiptf.T0_, addr);
    inletDir_.rmap(tiptf.inletDir_, addr);
}


void temperatureDirectedInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField& C = patch().Cf();
    vectorField rotationVelocity = omega_ ^ C;
    vectorField inletDirCCS = inletDir_/mag(inletDir_);
    vectorField inletDir(inletDirCCS);

    if (cylindricalCS_)
    {
        scalar radius;
        forAll(C, facei)
        {
            radius = sqrt(sqr(C[facei].x())+ sqr(C[facei].y()));
            inletDir[facei].z() = inletDirCCS[facei].z();

            if (radius > 0.0)
            {
                inletDir[facei].x() =
                    (C[facei].x()*inletDirCCS[facei].x()
                    - C[facei].y()*inletDirCCS[facei].y())/radius;

                inletDir[facei].y() =
                    (C[facei].y()*inletDirCCS[facei].x()
                    + C[facei].x()*inletDirCCS[facei].y())/radius;
            }
            else
            {
                inletDir[facei].x() = 0.0;
                inletDir[facei].y() = 0.0;
            }
        }
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        vectorField n = patch().nf();
        scalarField ndmagS = (n & inletDir_)*patch().magSf();

        operator==((inletDir_*phip/stabilise(ndmagS,SMALL))-rotationVelocity);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& Tp =
            patch().lookupPatchField<volScalarField, scalar>(TName_);

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        volScalarField Cp = thermo.Cp();

        const fvPatchField<scalar>& Cpp =
            patch().patchField<volScalarField, scalar>(Cp);

        operator==((inletDir*sqrt(2.0*Cpp*max(T0_-Tp,SMALL)))-rotationVelocity);
    }
    else
    {
        FatalErrorIn
        (
            "temperatureDirectedInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void temperatureDirectedInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (TName_ != "T")
    {
        os.writeKeyword("T") << TName_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("cylindricalCS") << cylindricalCS_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega")<< omega_ << token::END_STATEMENT << nl;
    T0_.writeEntry("T0", os);
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void temperatureDirectedInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(inletDir_*(inletDir_ & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    temperatureDirectedInletVelocityFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
