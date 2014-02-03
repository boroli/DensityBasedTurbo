    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

// Requires a Tpsi thermo class:
//    volScalarField& T = thermo.T();

    volScalarField& h = thermo.h();    
    volScalarField& p = thermo.p();       
    const volScalarField& T = thermo.T();
    
    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "compressibleCreatePhi.H"

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Create Riemann solver\n" << endl;
    godunovFlux<roeMultiDensityALEFlux> Godunov(p, U, T, thermo, turbulence());
    
    Info<< "Create kappa\n" << endl;
    volScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.Cp()/thermo.Cv()
    );
    
    Info<< "Create rhoU\n" << endl;
    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    Info<< "Create rhoE\n" << endl;
    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*(h + 0.5*magSqr(U))-p// + turbulence->k()) - p
    );

    Info<< "Create local time-step\n" << endl;    
    localTimeStep localTimeStep(U, thermo, turbulence());

    // relative velocity
    volScalarField Ur
    (
        IOobject
        (
            "Ur",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar("Ur", dimVelocity, 1.0)
    );
    // heat velocity
    volScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar("c", dimVelocity, 0.0)
    );    
        
    // pResidual
    volScalarField pResidual
    (
        IOobject
        (
            "pResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedScalar("pResidual", dimensionSet(1,-1,-3,0,0,0,0), scalar(1.0))
    );
    // UResidual
    volVectorField UResidual
    (
        IOobject
        (
            "UResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedVector("UResidual", dimensionSet(0,1,-2,0,0,0,0), vector::one)
    );
    // TResidual
    volScalarField TResidual
    (
        IOobject
        (
            "TResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar("TResidual", dimensionSet(0,0,-1,1,0,0,0), scalar(1.0))
    );    
    IOField<scalar> physDeltaT
    (
        IOobject
        (
            "physDeltaT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.0
    );
