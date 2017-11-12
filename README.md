# OFSolvers  
A project for various programs under OpenFOAM-5.0 framework. **WIP**

## Solvers:  
- *reactingDiffusionFoam*: A solver for reacting flows with multispecies diffusion using Hirshfelder and Curtiss approximation. Modification of *reactingFoam*;
- *ignitionFoam*: A solver with Energy Deposition model (see *Lacaze et al., 2009*) for spark/laser ignition modelling. Modification of *reactingFoam*

## Boundary conditions  
- *chokedInletPressure*: A pressure boundary condition for choked inlet flow. Permit to switch automatically between sonic (where pressure must be provided) and subsonic (where pressure is taken from internal field) regimes

## Examples

- *counterFlowFlame2D_ignitionFoam* - tutorial for *ignitionFoam* with an ignition spot modelling.

## Future work

* choked inlet boundary conditions for velocity and temperature: `chokedInletVelocity` and `chokedInletTemperature`

##Bibliography:  
* T. Poinsot and D. Veynante. *Theoretical and Numerical Combustion*, 3rd edition
* F. A. Williams. *Combustion Theory*
* G. Lacaze, B. Cuenot, T. Poinsot and M. Oschwald. *Large Eddy Simulation of laser ignition and compressible reacting flow in a rocket-like configuration* (2009)