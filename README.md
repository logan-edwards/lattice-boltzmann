# lattice-boltzmann
A simple implementation of the Lattice Boltzmann Method with solutions to simple problems in computational fluid dynamics (CFD).

lbm.c contains functions for each step of the actual LBM algorithm, while the demo files contain implementations used during development to verify behavior of the solver. Boundary conditions are handled on a per-problem basis between the collision step and the streaming step with custom handling depending on the geometry.
