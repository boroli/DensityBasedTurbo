//  creates Stress Mesh

    autoPtr<dynamicFvMesh> solidMeshPtr
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                "solid",
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    dynamicFvMesh& stressMesh = solidMeshPtr();
