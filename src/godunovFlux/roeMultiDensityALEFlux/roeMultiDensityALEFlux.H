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
    roeALEFlux

Description
    Implementation of inviscid numerical fluxes based on Roe & Pike FDS scheme
    with Harten's entropy fix, in an ALE formulation.
    Roe flux with preconditioning according to Weiss & Smith:
    AIAA Journal, Vol. 33, No. 11, 1995

    ALE formulations:
    AIAA-98-2411, Darracq et.al
    AIAA-2006-2825, Gönc et.al

Author
    Hrvoje Jasak All rights reserved.
    Aleksandar Jemcov  All rights reserved.
    Oliver Borm  All rights reserved.
    Sebastian Saegeler  All rights reserved.

SourceFiles
    roeMultiDensityALEFlux.H

\*---------------------------------------------------------------------------*/

#ifndef roeMultiDensityALEFlux_H
#define roeMultiDensityALEFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class roeMultiDensityALEFlux Declaration
\*---------------------------------------------------------------------------*/

class roeMultiDensityALEFlux
{

public:

    // Constructors

        //- Construct from components
        roeMultiDensityALEFlux(){};

    // Destructor - default

    // Member Functions

        //- Update fluxes based on current left and right states
        void evaluateFlux
        (
            scalar& rhoFlux,
            vector& rhoUFlux,
            scalar& rhoEFlux,
            const scalar pLeft,
            const scalar pRight,
            const vector ULeft,
            const vector URight,
// //
//     const scalar TLeft,
//     const scalar TRight,
// //
            const scalar rhoLeft,
            const scalar rhoRight,
// //
            const scalar visVelLeft,
            const scalar visVelRight,
// //
     const scalar RLeft,
     const scalar RRight,
//     const scalar CvLeft,
//     const scalar CvRight,
// //
            const scalar kappaLeft,
            const scalar kappaRight,
// //
            const vector Sf,
            const scalar magSf,
            const vector dotX,
            const scalar Konstant
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
