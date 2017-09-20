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
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);
    volScalarField Dki(0.0*Dk[0]);

    // Mixture molecular weight
    Wm = composition.W();
    // Dumping correction velocity before the new iteration
    Vc = dimensionedVector("Vc",dimVelocity,vector(0,0,0));
    // Compute correction velocity as a sum of all species
    forAll(Y, i)
    {
        scalar Wk = composition.W(i);
        // Define molar fraction for specie Yi
        volScalarField Xk = Y[i]*Wm/Wk;
        // Compute diffusion coefficient for species i
        // Loop sorting because of different order of species in Y[] and species[]
        forAll(species, j)
        {
            if (species[j] == Y[i].name())
            {
                Dk[j] = thermo.alpha()/rho/Le[j];
                Dki = Dk[j];
            }
        }

        Vc += Dki*Wk/Wm*fvc::grad(Xk);
    }
    // Compute correction velocity flux
    surfaceScalarField phiVc = linearInterpolate(rho*Vc) & mesh.Sf();
    // Divergence scheme for correction velocity
    tmp<fv::convectionScheme<scalar>> mvConvectionCorr
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phiVc,
            mesh.divScheme("div(phiVc,Yi_h)")
        )
    );

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            forAll(species, j)
            {
                if (species[j] == Yi.name())
                {
                    Dk[j] = thermo.alpha()/rho/Le[j];
                    Dki = Dk[j];
                }
            }
            volScalarField expDiffTerm = fvc::laplacian(rho*Dki*Yi/Wm,Wm);

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
            // Additional convection term due to correction velocity
              + mvConvectionCorr->fvmDiv(phiVc, Yi)
            /* Deprecated term because Fick's law and Sh=1 is not valid anymore
              - fvm::laplacian(turbulence->muEff(), Yi)*/
            // Adding instead turbulent flux based on mut() and turbulent Schmidt number
              - fvm::laplacian(turbulence->mut()/Sht,Yi)
             ==
                reaction->R(Yi)
            /* Additional terms due to multispecies diffusion
            Decoupled in order to construct discretization around Yi */
              + fvm::laplacian(rho*Dki, Yi)
              + expDiffTerm
            //------------------------------------------------//
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
