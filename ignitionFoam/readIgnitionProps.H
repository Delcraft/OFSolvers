IOdictionary ignitionProperties
(
    IOobject
	(
    	    "ignitionProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
     )
);

const wordList spotNames(ignitionProperties.toc());
 
    PtrList<dimensionedScalar> eps_(spotNames.size());
    PtrList<dimensionedScalar> sigmas_(spotNames.size());
    PtrList<dimensionedScalar> sigmat_(spotNames.size());
    PtrList<dimensionedScalar> tmax_(spotNames.size());
    PtrList<dimensionedVector> location_(spotNames.size());
 
    forAll(spotNames, i)
    {
        const word& spotName = spotNames[i];
        const dictionary& subDict = ignitionProperties.subDict(spotName);
 
        eps_.set
        (
             i,
             new dimensionedScalar(subDict.lookup("eps"))
        );
        sigmas_.set
        (
             i,
             new dimensionedScalar(subDict.lookup("sigmas"))
        );
        sigmat_.set
        (
             i,
             new dimensionedScalar(subDict.lookup("sigmat"))
        );
        tmax_.set
        (
             i,
             new dimensionedScalar(subDict.lookup("tmax"))
        );
        location_.set
        (
             i,
             new dimensionedVector(subDict.lookup("location"))
        );
        Info<< "eps = " << eps_[i] << endl;
        Info<< "sigmas = " << sigmas_[i] << endl;
        Info<< "sigmat = " << sigmat_[i] << endl;
        Info<< "tmax = " << tmax_[i] << endl;
        Info<< "location = " << location_[i] << endl;
    }
