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

#include "dynamicSolidFvMesh.H"
#include "addToRunTimeSelectionTable.H"
// #include "volFields.H"
// #include "mathematicalConstants.H"

// #include "pointMesh.H"
// #include "pointFields.H"

#include "volPointInterpolation.H"
#include "ZoneIDs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicSolidFvMesh, 0);

addToRunTimeSelectionTable(dynamicFvMesh, dynamicSolidFvMesh, IOobject);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicSolidFvMesh::dynamicSolidFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    fvMesh_(refCast<const fvMesh>(*this)),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement",
            io.time().timeName(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(fvMesh_),
        dimensionedVector
        (
            "pointDisplacement",
            pointDisplacement_.dimensions(),
            vector::zero
        )
    ),
    stationaryPoints_
    (
            IOobject
            (
                "points",
                io.time().constant(),
                meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
    ),
//     pointLocation_(NULL),
    points0_
    (
        pointIOField
        (
            IOobject
            (
                "points",
                io.time().constant(),
                meshSubDir,
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        )
    )
{

// //     if (io.headerOk())
// //     {
//         pointLocation_.reset
//         (
//             new pointVectorField
//             (
//                 io,
//                 pointMesh::New(fvMesh_)
//             )
//         );
// 
//         if (debug)
//         {
//             Info<< "displacementLaplacianFvMotionSolver :"
//                 << " Read pointVectorField "
//                 << io.name() << " to be used for boundary conditions on points."
//                 << nl
//                 << "Boundary conditions:"
//                 << pointLocation_().boundaryField().types() << endl;
//         }
// //     }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynamicSolidFvMesh::~dynamicSolidFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool dynamicSolidFvMesh::update()
{
//     // interpolates volume fields to point fields for moving mesh
    const volPointInterpolation& cpi = volPointInterpolation::New(*this);

//     const objectRegistry& registry = *this;

    const volVectorField& Usolid =
        lookupObject<volVectorField>("U");

    // takes the data on the solid patch (cell centres)
    // and interpolates it into points
    pointVectorField solidPointsDispl =
//         cpi.interpolate(Usolid - Usolid.oldTime());
        cpi.interpolate(Usolid);
	
// 	pointVectorField solidPointsDispl =
// 		volPointInterpolation::New(*this).interpolate(Usolid);

//     volPointInterpolation::New(*this).interpolate
//     (
//         Usolid,
//         pointDisplacement_
//     );

// 	volPointInterpolation::New(*this).interpolate
// 	(
// 		Usolid,
// 		pointDisplacement_
// 	);
	
// 	const vectorField& sldDispl = solidPointsDispl.internalField();

// 	pointVectorField movingPoints(allPoints().size(),vector::zero);
	
	// TODO
// 	pointDisplacement_.internalField() =
// 		stationaryPoints_
// 		+ solidPointsDispl.internalField();
// 	
// 	pointDisplacement_.correctBoundaryConditions();

//         pointLocation_().internalField() =
//             points0_
//           + pointDisplacement_.internalField();
// 
//         pointLocation_().correctBoundaryConditions();

//     // Retrieve the cell zone Names
//     const wordList pointZoneNames = pointZones().names();
//     
//     Info << pointZoneNames << pointZoneNames << endl;
// 
//     // Set the points
// //     movingPointsPtr_ = new vectorField(allPoints().size(),vector::zero);
// 
//     forAll(pointZoneNames,cellZoneI)
//     {
// 	if (pointZoneNames[cellZoneI] == "allMyPoints")
// 	{
// 		const labelList& pointAddr =
// 			pointZones()[pointZones().findZoneID(pointZoneNames[cellZoneI])];
// 			
// 		Info << "pointAddr.size(): " << pointAddr.size() << endl;
// 		
// 
// 	}
//     }
// 		forAll(movingPoints,pointI)
// 		{
// 			movingPoints[pointI] =
// 				stationaryPoints_[pointI]
// 				+ sldDispl[pointI];
// 			
// 	// 		reduce(movingPoints[pointI], sumOp<vector>());
// 		}

//  movingPoints[curFace[pointI]] =

// 	solidPointsDispl.correctBoundaryConditions();
	
// 	Info << "stationaryPoints_.size(): " << stationaryPoints_.size() << endl;


// 	points0_.transfer(solidPointsDispl.internalField());

//         tmp<pointField> tcurPoints
//         (
//             points0_ + pointDisplacement_.internalField()
//         );


    fvMesh::movePoints(stationaryPoints_);

    // Mesh motion only - return false
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
