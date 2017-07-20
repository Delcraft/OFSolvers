IOdictionary speciesDict
(
    IOobject
    (
        "transportProperties",    // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
        IOobject::MUST_READ_IF_MODIFIED,    // must exist, otherwise failure
        IOobject::NO_WRITE      // dict is only read by the solver
    )
);

//const wordList specieNames(speciesDict.toc());

const wordList species(speciesDict.lookup("species"));

const dictionary& LeDict(speciesDict.subDict("HirshfelderCurtiss"));

List<dimensionedScalar> Le(species.size());

//List<dimensionedScalar> Le_(specieNames.size());

Info << "Reading mass diffusion settings for all species" << endl;

//Info << species << endl;

forAll(species, i)
{
    word thisSpecie = species[i];

    if (LeDict.found(species[i]))
    {
        Le[i] = dimensionedScalar(LeDict.lookup(thisSpecie));
    }
    else
    {
        FatalErrorIn
        (
            "const dictionary& LeDict(speciesDict.subDict(""HirshfelderCurtiss""))"
        ) << "Specie " << species[i] << " not found" << endl
        << "Species used in this model are: " << species << endl
        << exit(FatalError);
    }
//    Le = subDict.lookup("Le");

    Info << "Le number for specie " << thisSpecie << ": " << Le[i].value() << endl;
}

List<autoPtr<volScalarField>> Dk(species.size());
forAll(species, i)
{
    word nameDk ("Dk_"+species[i]);

    forAll(Le, j)
    {
        if (Le[j].name() == species[i])
        {
            Dk[i].reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        nameDk,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    thermo.alpha() / rho / Le[j].value()
                )
            );
        }
    }
}
/*
forAll(species, i)
{
    forAll(Le, j)
    {
        if (Le[j].name() == species[i])
        {
            Info << "Test output: " << Le[j] << endl;
        }
    }
}
*/