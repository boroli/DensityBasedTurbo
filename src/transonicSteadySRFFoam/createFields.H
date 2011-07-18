    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

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
    godunovFlux<hllcALEFlux> Godunov(p, U, T, thermo, turbulence());
//     godunovFlux<AUSMplusALEFlux> Godunov(p, U, T, thermo, turbulence());
//     godunovFlux<roeALEFlux> Godunov(p, U, T, thermo, turbulence());
//     godunovFlux<rusanovALEFlux> Godunov(p, U, T, thermo, turbulence());   

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
        rho*(h + 0.5*magSqr(U) + turbulence->k()) - p
    );

    localTimeStep localTimeStep(U, thermo, turbulence());

    Info<< "\nCreating SRF model\n" << endl;
    autoPtr<SRF::SRFModel> SRF
    (
        SRF::SRFModel::New(U)
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
        1
    );
