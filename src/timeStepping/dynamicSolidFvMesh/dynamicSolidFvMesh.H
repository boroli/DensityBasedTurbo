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
    dynamicSolidFvMesh

Description
    Automatic simplified mesh motion for "box-in-mesh" cases.  Here,
    a direction of motion is defined, together with motion amplitude
    and frequency.  The domain is separated into three parts, where
    the middle part moves according to the prescribed motion law.
    Parts of the mesh before and after the obstacle are scaled.

SourceFiles
    dynamicSolidFvMesh.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef dynamicSolidFvMesh_H
#define dynamicSolidFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicSolidFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicSolidFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- The fvMesh to be moved
        const fvMesh& fvMesh_;

        //- Point motion field
        mutable pointVectorField pointDisplacement_;

        //- Stationary point location
        pointIOField stationaryPoints_;

        //- Optionally read point-position field. Used only for position
        //  boundary conditions.
//         mutable autoPtr<pointVectorField> pointLocation_;

        //- Reference point field
        pointField points0_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicSolidFvMesh(const dynamicSolidFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicSolidFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicSolidFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicSolidFvMesh(const IOobject& io);


    // Destructor

        virtual ~dynamicSolidFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
