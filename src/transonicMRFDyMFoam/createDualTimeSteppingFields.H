    label noScalarFields = 2;
    label noVectorFields = 1;

    wordList scalarFieldNames;
    wordList vectorFieldNames;

    PtrList<volScalarField> scalarFieldList;
    PtrList<volVectorField> vectorFieldList;

    if(numberSubCycles > 1)
    {
        Info<< "Creating oldTime fields for dual-time stepping\n" << endl;

        scalarFieldNames.setSize(noScalarFields);
        vectorFieldNames.setSize(noVectorFields);

        wordList volScalarFields = mesh.names("volScalarField");
        wordList volVectorFields = mesh.names("volVectorField");

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the turbulenceModel is created otherwise the dictionary is
    // entered in the database twice
    {
        IOdictionary turbulenceCoeffDict
        (
            IOobject
            (
                "turbulenceProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if(turbulenceCoeffDict.found("turbulentFields"))
        {
            ITstream turbulentFieldsInput
                = turbulenceCoeffDict.lookup("turbulentFields");

            word turbulentField;

            forAll(turbulentFieldsInput, fieldi)
            {
                turbulentFieldsInput >> turbulentField;

                forAll(volScalarFields, scalarFieldi)
                {
                    word scalField = volScalarFields[scalarFieldi];
                    if(scalField == turbulentField)
                    {
                        scalarFieldNames.setSize(scalarFieldNames.size()+1);
                        scalarFieldNames[scalarFieldNames.size()-1]
                            = turbulentField;
                    }
                }
                forAll(volVectorFields, vectorFieldi)
                {
                    word vecField = volVectorFields[vectorFieldi];
                    if(vecField == turbulentField)
                    {
                        vectorFieldNames.setSize(vectorFieldNames.size()+1);
                        vectorFieldNames[vectorFieldNames.size()-1]
                            = turbulentField;
                    }
                }
            }
        }
    }

        scalarFieldNames[0] = rho.name();
        scalarFieldNames[1] = rhoE.name();

        vectorFieldNames[0] = rhoU.name();

//     Info << "scalarFieldNames " << scalarFieldNames << endl;
//     Info << "vectorFieldNames " << vectorFieldNames << endl;

        scalarFieldList.setSize(2*scalarFieldNames.size());
        vectorFieldList.setSize(2*vectorFieldNames.size());

        // TODO: check second order restart in physical time
        word oldName = "Old";

        // you have to make sure, you are providing the correct values and time
        // steps!!!!
        forAll(scalarFieldNames, fieldi)
        {
            const volScalarField& tmpField =
                mesh.lookupObject<volScalarField>(scalarFieldNames[fieldi]);

            scalarFieldList.set
            (
                fieldi,
                new volScalarField
                (
                    scalarFieldNames[fieldi]+oldName,
                    tmpField
                )
            );

            scalarFieldList.set
            (
                fieldi+scalarFieldNames.size(),
                new volScalarField
                (
                    IOobject
                    (
                        scalarFieldNames[fieldi]+oldName+oldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    tmpField.oldTime()
                )
            );
        }

        forAll(vectorFieldNames, fieldi)
        {
            const volVectorField& tmpField =
                mesh.lookupObject<volVectorField>(vectorFieldNames[fieldi]);

            vectorFieldList.set
            (
                fieldi,
                new volVectorField
                (
                     vectorFieldNames[fieldi]+oldName,
                     tmpField
                )
            );

            vectorFieldList.set
            (
                fieldi+vectorFieldNames.size(),
                new volVectorField
                (
                    IOobject
                    (
                        vectorFieldNames[fieldi]+oldName+oldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    tmpField.oldTime()
                )
            );
        }
    }
