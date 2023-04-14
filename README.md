# Analysis of rescattering effects in $3\pi$ final states

Decays into three particles are often described in terms of two-body resonances and a non-interacting spectator particle.  To go beyond this simplest isobar model, crossed-channel rescattering effects need to be accounted for.
We quantify the importance of these rescattering effects in three-pion systems for different decay masses and angular-momentum quantum numbers.
We provide the amplitude decompositions for four decay processes with total $J^{PC} = 0^{--}$, $1^{--}$, $1^{-+}$, and $2^{++}$, 
all of which decay predominantly as $\rho\pi$ states.
Two-pion rescattering is described in terms of an Omn{\`e}s function, which incorporates the $\rho$ resonance. 
Inclusion of crossed-channel effects is achieved by solving the Khuri--Treiman integral equations. 
The unbinned log-likelihood estimator is used to determine the significance of the rescattering effects beyond two-body resonances;
we compute the minimum number of events necessary to unambiguously find these in future Dalitz-plot analyses.
Kinematic effects that enhance or dilute the rescattering are identified for 
the selected set of quantum numbers and various masses.

## Code

### Additional exploration (Julia)

- [`N-001-cartesian-pions.jl`](crosscheck/N-001-cartesian-pions.jl) Relation of the physical decay to the isospin amplitudes in cartesian coordinates
- [`N-002-covariant3pi.jl`](crosscheck/N-002-covariant3pi.jl) Evaluation of the covariant structures in the s-channel
- [`N-003-omega.jl`](crosscheck/N-003-omega.jl) early-stage MC experiments on discriminating models
