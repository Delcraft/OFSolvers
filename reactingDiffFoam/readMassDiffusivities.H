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

const wordList species(speciesDict.lookup("species"));

const dictionary& LeDict(speciesDict.subDict("HirshfelderCurtiss"));

List<dimensionedScalar> Le(species.size());

Info << "Using Hirshfelder and Curtiss approximation"
     << " for multispecies diffusion flux" << endl;

Info << "Reading mass diffusion settings for all species" << endl;

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

    Info << "Le number for specie " << thisSpecie << ": " << Le[i].value() << endl;
}

scalar Sht = 1.0;

if (LeDict.found("Sht"))
{
    Sht = readScalar(LeDict.lookup("Sht"));
    Info << "Sht = " << Sht << endl;
}
else
{
    Info << "Sht is not found in " << speciesDict.name() << ". Taking Sht = 1.0" << endl;
}

PtrList<volScalarField> Dk(species.size());
forAll(species, i)
{
    word nameDk ("Dk"+species[i]);

    Dk.set
    (
        i,
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
//                    mesh,
//                    dimensionedScalar(nameDk,dimensionSet(0,2,-1,0,0,0,0),0.0)
            thermo.alpha() / rho / Le[i]
        )
    );
}

// Correction velocity
volVectorField Vc
(
    IOobject
    (
        "Vc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Vc",dimVelocity,vector(0,0,0))
);