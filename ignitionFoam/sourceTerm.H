const scalar& PI = constant::mathematical::pi;

dimensionedScalar t
(
"t",
dimensionSet(0,0,1,0,0,0,0),
runTime.value()
);
//Info << "current time for ST: " << t << endl;

//Setting all values of ST from previous time to zero
ST = dimensionedScalar("init_zero",dimensionSet(1,-1,-3,0,0,0,0),0.0);

forAll(spotNames,i)
{
        dist = mag(mesh.C() - location_[i]);
        ST += eps_[i]/pow(sigmas_[i],3)/sigmat_[i]/4/sqr(PI)
        *exp(-0.5*sqr(dist/sigmas_[i]))*exp(-0.5*sqr((t - tmax_[i])/sigmat_[i]));
}

