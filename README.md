# numericalSolvers
4th-order Runge Kutta &amp; Explicit Euler for the generalisation of a mass spring system.

parameters.txt holds the system parameters and initial state quantities of the masses.
The first line of parameters.txt looks like **a b c d e** where:y
**a** is the binary value for the time integration scheme (0 for Explicit Euler, 1 for 4th order Runge-Kutta)
**b** is the upper time bound on which the integration scheme acts 
**c** is the time step
**d** is the spring constant of the system
**e** is the damping factor applied to each mass individually (i.e. it only depends on the velocity of the mass it acts on)

All other lines hold 3 values, **f g h** where:
**f** is the mass of the i'th line
**g** is the initial displacement of the i'th line and,
**h** is the initial velocity of the i'th mass.

The parameters.txt file MUST be laid out as above, otherwise the program will not function correctly.
The code will terminate early if any of the following are negative: damping, stiffness, timestep, upper time bound or mass.
If **a** is not recognised, the code will also terminate early.

%TODO - can you include instructions for complimation?