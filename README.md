# Routines for the computations of periodic orbits using Melnikov and Newton method in a piecewise-smooth system (the rocking block)
Collection of C and bash routines for the computation of periodic orbits via Melnikov method in a piecewise-smooth system
This software is the one used to compute the numerical results presented in 
<a href="http://dx.doi.org/10.1137/110850359">
The Melnikov method and subharmonic orbits in a piecewise-smooth system</a> .

**List of files, routines and their goals**
* system.dat contains system parameters: n m y0 t0 delta rho omega, where
<ul>
<li> nT is the period of the desired periodic orbit</li>
<li> 2*m the number of impacts</li>
<li>y0 is the initial seed for the Newton method, such that (0,y0) is the initial condition for an nT period orbit of the unperturbed system </li>
<li> t0 is the zero of the modified Melnikov function </li>
<li>delta is perturbation parameter </li>
<li> rho is the ratio between the forcing and the restitution coefficient </li>
<li> omega is the frequency of the forcing </li>
</ul>
* find_nm_ini_cond computes the initial conditions, (x0,t0) for a periodic orbit whose parameters are given in system.dat
* solve_wt0 computes a zero of the Melnikov function
* find_delta_max.sh and follow_orbit.sh perform a continuation method by increasing delta until the periodic orbit bifurcates
* find_existence_regions.sh computes the existence region for an n,m-periodic orbit in the epsilon-r parameter space
