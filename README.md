# Analysis of rescattering effects in $3\pi$ final states

[![GitHub page](https://img.shields.io/badge/GitHub-README.md-yellowgreen)](https://github.com/mmikhasenko/rescattering-3pi-2023-001)
[![arXiv article](https://img.shields.io/badge/article-%20hep--ph%3A2212.11767-brightgreen)](https://inspirehep.net/literature/2617378)


Decays into three particles are often described in terms of two-body resonances and a non-interacting spectator particle.
To go beyond this simplest isobar model, crossed-channel rescattering effects need to be accounted for.
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

The main numerics is done with c++/python code which can be found in the repository [`HISKP-ph/khuri_treiman_solver`](https://github.com/HISKP-ph/khuri_treiman_solver).

### Julia code

- [`KTMC.jl/scr`](https://github.com/mmikhasenko/rescattering-3pi-2023-001/tree/main/KTMC.jl/src) source code the Khuri-Treiman Monte-Carlo project
- [`KTMC.jl/scripts`](https://github.com/mmikhasenko/rescattering-3pi-2023-001/tree/main/KTMC.jl/scripts) several exploratory scripts for the project
- [`KTMC.jl/plots`](https://github.com/mmikhasenko/rescattering-3pi-2023-001/tree/main/KTMC.jl/plots) preliminary plots

### Pluto notebooks (Julia)

- [`N-001-cartesian-pions`](crosscheck/N-001-cartesian-pions.html) [[code]](https://github.com/mmikhasenko/rescattering-3pi-2023-001/blob/main/crosscheck/N-001-cartesian-pions.jl) Relation of the physical decay to the isospin amplitudes in cartesian coordinates
- [`N-002-covariant3pi`](crosscheck/N-002-covariant3pi.html) [[code]](https://github.com/mmikhasenko/rescattering-3pi-2023-001/blob/main/crosscheck/N-002-covariant3pi.jl) Evaluation of the covariant structures in the s-channel
- [`N-003-omega`](crosscheck/N-003-omega.html) [[code]](https://github.com/mmikhasenko/rescattering-3pi-2023-001/blob/main/crosscheck/N-003-omega.jl) early-stage MC experiments on discriminating models
