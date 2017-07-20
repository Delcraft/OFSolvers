tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    dQ = reaction->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {

        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
/*
            forAll(Dk, j)
            {
                if (Dk[j].value() == Yi.name())
                {
                    Info << "Dk name for specie " << Yi.name() << ": " << Dk[j].value() << endl;
                }
            }
*/
//            Info << "Le number for specie " << Yi.name() << ": " << Le(Yi.name()) << endl;

//            volScalarField expDiffTerm = fvc::laplacian(rho*Dk*Yi/Wm,Wm);

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
