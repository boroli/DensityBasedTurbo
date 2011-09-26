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
    Foam::BarthJespersenSlopeMultiLimiter

Description
    Class with limiter function which returns the limiter
    BarthJespersenSlopeMulti

Author
    Oliver Borm  All rights reserved.

SourceFiles
    BarthJespersenSlopeMulti.H

\*---------------------------------------------------------------------------*/

#ifndef BarthJespersenSlopeMulti_H
#define BarthJespersenSlopeMulti_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BarthJespersenSlopeMultiLimiter Declaration
\*---------------------------------------------------------------------------*/

class BarthJespersenSlopeMultiLimiter
{

public:

    BarthJespersenSlopeMultiLimiter(Istream&)
    {}

    // Barth-Jespersen Slope Limiter
    // Implementation acccording to Blazek Chap.5 / p. 166
    inline scalar limiter
    (
        const scalar cellVolume,
        const scalar deltaOneMax,
        const scalar deltaOneMin,
        scalar deltaTwo
    )
    {
        deltaTwo = stabilise(deltaTwo,VSMALL);
        return min(max(max(deltaOneMax/deltaTwo,0),max(deltaOneMin/deltaTwo,0)),1);
    }

    inline scalar limiter
    (
        const scalar cellVolume,
        const scalar faceFlux,
        const scalar deltaOneMax,
        const scalar deltaOneMin,
        const vector gradPhiP,
        const vector gradPhiN,
        const vector& d
    )
    {
        return limiter(cellVolume,deltaOneMax,deltaOneMin,gradPhiP & d);
    }

    inline vector limiter
    (
        const scalar cellVolume,
        const scalar faceFlux,
        const vector deltaOneMax,
        const vector deltaOneMin,
        const tensor gradPhiP,
        const tensor gradPhiN,
        const vector& d
    )
    {
        vector deltaTwo = d & gradPhiP;
        return vector
        (
            limiter(cellVolume,deltaOneMax[0],deltaOneMin[0],deltaTwo[0]),
            limiter(cellVolume,deltaOneMax[1],deltaOneMin[1],deltaTwo[1]),
            limiter(cellVolume,deltaOneMax[2],deltaOneMin[2],deltaTwo[2])
        );
    }
/*
//     inline vector limiter
    inline scalar limiter
    (
        const scalar cellVolume,
        const scalar faceFlux,
        const vector deltaOneMax,
        const vector deltaOneMin,
        const tensor gradPhiP,
        const tensor gradPhiN,
        const vector& d
    )
    {
        vector deltaTwo = gradPhiP & d;
//         scalar myLimiter = min(min
//         (
//             limiter(cellVolume,deltaOneMax[0],deltaOneMin[0],deltaTwo[0]),
//             limiter(cellVolume,deltaOneMax[1],deltaOneMin[1],deltaTwo[1])),
//             limiter(cellVolume,deltaOneMax[2],deltaOneMin[2],deltaTwo[2])
//         );
//         return myLimiter;

        return mag(vector
        (
            limiter(cellVolume,deltaOneMax[0],deltaOneMin[0],deltaTwo[0]),
            limiter(cellVolume,deltaOneMax[1],deltaOneMin[1],deltaTwo[1]),
            limiter(cellVolume,deltaOneMax[2],deltaOneMin[2],deltaTwo[2])
        ))/3.0;

//         return vector(myLimiter,myLimiter,myLimiter);

    }
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
