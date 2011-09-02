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

\*---------------------------------------------------------------------------*/

#include "outletInletBackflowFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
outletInletBackflowFvPatchField<Type>::outletInletBackflowFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF)
{
    this->refValue() = *this;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
outletInletBackflowFvPatchField<Type>::outletInletBackflowFvPatchField
(
    const outletInletBackflowFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
outletInletBackflowFvPatchField<Type>::outletInletBackflowFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF)
{
    this->refValue() = Field<Type>("outletValue", dict, p.size());
/*
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {*/
        fvPatchField<Type>::operator=(this->refValue());
//     }

    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
outletInletBackflowFvPatchField<Type>::outletInletBackflowFvPatchField
(
    const outletInletBackflowFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf)
{}


template<class Type>
outletInletBackflowFvPatchField<Type>::outletInletBackflowFvPatchField
(
    const outletInletBackflowFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void outletInletBackflowFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

//     const fvsPatchField<scalar>& phip = this->patch().lookupPatchField
//     (
//         "phi",
//         reinterpret_cast<const surfaceScalarField*>(0),
//         reinterpret_cast<const scalar*>(0)
//     );
//
//     this->valueFraction() = pos(phip - SMALL);

    const fvPatchField<vector>& Up = this->patch().lookupPatchField
    (
        "U",
        reinterpret_cast<const volVectorField*>(0),
        reinterpret_cast<const vector*>(0)
    );

    this->valueFraction() = pos((this->patch().nf() & Up) - SMALL);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void outletInletBackflowFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->refValue().writeEntry("outletValue", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
