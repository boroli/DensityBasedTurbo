    if (numberSubCycles > 1)
    {
        // you have to make sure, you are providing the correct values and time
        // steps!!!!
        forAll(scalarFieldNames, fieldi)
        {
            const volScalarField& tmpField =
                mesh.lookupObject<volScalarField>(scalarFieldNames[fieldi]);

            scalarFieldList[fieldi+scalarFieldNames.size()] =
                scalarFieldList[fieldi];
            scalarFieldList[fieldi+scalarFieldNames.size()].timeIndex() =
                scalarFieldList[fieldi].timeIndex();

            scalarFieldList[fieldi] = tmpField;
            scalarFieldList[fieldi].timeIndex() = tmpField.time().timeIndex();
        }

        // you have to make sure, you are providing the correct values and time
        // steps!!!!
        forAll(vectorFieldNames, fieldi)
        {
            const volVectorField& tmpField =
                mesh.lookupObject<volVectorField>(vectorFieldNames[fieldi]);

            // use copy constructor, as this saves the correct timeIndex
            vectorFieldList[fieldi+vectorFieldNames.size()] =
                vectorFieldList[fieldi];
            vectorFieldList[fieldi+vectorFieldNames.size()].timeIndex() =
                vectorFieldList[fieldi].timeIndex();

            vectorFieldList[fieldi] = tmpField;
            vectorFieldList[fieldi].timeIndex() = tmpField.time().timeIndex();
        }
    }
