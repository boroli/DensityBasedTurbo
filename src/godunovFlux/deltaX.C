// calculation the caracteristic length in each cell
// necessary for viscous velocity
// (has to be consistent with local-TS method)

    volScalarField deltaX
    (
        IOobject
        (
            "deltaXGodunov",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("deltaXGodunov", dimLength, HUGE),
        zeroGradientFvPatchScalarField::typeName
    );


    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector normal = Sf[faceI]/magSf[faceI]; //face normal vector

        deltaX[own] = 
            min(deltaX[own], mag((faceCenter[faceI]-cellCenter[own]) & normal));

        deltaX[nei] =
            min(deltaX[nei], mag((faceCenter[faceI]-cellCenter[nei]) & normal));
    }

    forAll(deltaX.boundaryField(), patchI)
    {

        const fvPatchScalarField& pp = deltaX.boundaryField()[patchI];
        const vectorField pFaceCenter = pp.patch().Cf();

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];

        const unallocLabelList& faceCells =  pp.patch().faceCells();

        forAll(pp, facei)
        {
            vector normal = pSf[facei]/pMagSf[facei];

            label own = faceCells[facei];

            deltaX[own] = min
            (
                deltaX[own],
                mag((pFaceCenter[facei]-cellCenter[own]) & normal)
            );
        }
    }
    
    deltaX.correctBoundaryConditions();  
