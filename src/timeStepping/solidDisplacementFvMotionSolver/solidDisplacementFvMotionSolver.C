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

#include "solidDisplacementFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidDisplacementFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        solidDisplacementFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidDisplacementFvMotionSolver::solidDisplacementFvMotionSolver
(
    const polyMesh& mesh,
    Istream& is
)
:
    displacementFvMotionSolver(mesh, is),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(fvMesh_),
        dimensionedVector
        (
            "pointDisplacement",
            dimLength,
            vector::zero
        )
    ),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellDisplacement",
            dimLength,
            vector::zero
        ) //,
//         cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(NULL),
//     diffusivityPtr_
//     (
//         motionDiffusivity::New(*this, lookup("diffusivity"))
//     ),
    frozenPointsZone_
    (
        found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID(lookup("frozenPointsZone"))
      : -1
    )
{
    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "solidDisplacementFvMotionSolver:" << nl
//             << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }


    if (io.headerOk())
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "solidDisplacementFvMotionSolver :"
                << " Read pointVectorField "
                << io.name() << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidDisplacementFvMotionSolver::
~solidDisplacementFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::solidDisplacementFvMotionSolver::curPoints() const
{
    volPointInterpolation::New(fvMesh_).interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    if (pointLocation_.valid())
    {
        if (debug)
        {
            Info<< "solidDisplacementFvMotionSolver : applying "
                << " boundary conditions on " << pointLocation_().name()
                << " to new point location."
                << endl;
        }

        pointLocation_().internalField() =
            points0()
          + pointDisplacement_.internalField();

        pointLocation_().correctBoundaryConditions();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                pointLocation_()[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(pointLocation_().internalField());

        return tmp<pointField>(pointLocation_().internalField());
    }
    else
    {
        tmp<pointField> tcurPoints
        (
            points0() + pointDisplacement_.internalField()
        );

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                tcurPoints()[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(tcurPoints());

        return tcurPoints;
    }
}


void Foam::solidDisplacementFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.points());

//     diffusivityPtr_->correct();
    pointDisplacement_.boundaryField().updateCoeffs();

    const objectRegistry& registry = fvMesh_;

    const volVectorField& Usolid =
        registry.lookupObject<volVectorField>("U");

//     cellDisplacement_ = Usolid;
    
    cellDisplacement_ = (Usolid - Usolid.oldTime());

//     // takes the data on the solid patch (cell centres)
//     // and interpolates it into points
//     pointVectorField solidPointsDispl =
// //         cpi.interpolate(Usolid - Usolid.oldTime());
//         cpi.interpolate(Usolid);

//     Foam::solve
//     (
//         fvm::laplacian
//         (
//             diffusivityPtr_->operator()(),
//             cellDisplacement_,
//             "laplacian(diffusivity,cellDisplacement)"
//         )
//     );
}


void Foam::solidDisplacementFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementFvMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
//     diffusivityPtr_.reset(NULL);
//     diffusivityPtr_ = motionDiffusivity::New(*this, lookup("diffusivity"));
}


// ************************************************************************* //
