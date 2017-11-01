# dilu_cham_modl
Code for simulating the effect of dilution on particle mass

Inputs
Environmental Variables
Temperature, pressure, time to solve over, molecular weight dry air
Bulk Particle Properties
Radius, number concentration
VBS Properties
Molecular weight, density, saturation vapour pressure, gas-phase concentration, mass fractions in particle phase, accommodation coefficient

Pre-time loop calculations
Gas-phase diffusion coefficient
The diffusion coefficient is a function of temperature, pressure, molecular weight of the partitioning component and air, and particle size.  It needs a correction applied when the aerosol is not in the continuum regime.
Corrected concentration above particle surface
The saturation concentration is found and then it is corrected for the curvature (Kelvin) effect and solute (Raoult) effect.


Time Loop
Estimate of concentration change based on difference between bulk vapour-phase concentration and concentration above particle surface, following Ficks first law (eq. 3 of Zaveri et al. 2008).  This is the governing equation for partitioning.

The code includes an option to use the partitioning equation with the effect of partitioning on particle radius through the mass transfer rate and the curvature effect accounted for.  There is also an option where the effect of partitioning on radius is presumed to be relatively small and so mass transfer rate and curvature effects are assumed constant.

Output
Mass fraction remaining in particle phase
