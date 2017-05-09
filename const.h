#ifndef CONST_H
#define CONST_H

/* growth rate
 */

#define growth_rate  0.0005
#define DISC_RADIUS  60
#define LENGTH_FIRST_LATTICE 8	/* Nb of cell in one row, has to be even! */
#define NU 0.001 /* 0.005 Differentiation prefactor */
#define MU 30	/* 6 Area factor */
#define MUP 1
#define MUT 50
#define RHO 1.0
#define RHOEXP 1.5
#define RHOFCT 2.0
#define DIVISION_RATE 0.0015 /* Define the probability of getting one cell division within one step */
#define PERIMETER_ON 1
#define BENDORINERT_TERM 0 /* 1 if there is a bending or inertia term in the hamiltonian */
#define AREA 1. /* default cell area */

/* Cutoff used to know whether the lattice is relaxed */
#define LEN_CUTOFF 1.0E-4
#define MAX_ITER 10000

#define WIN_SIZE 800, 800	/* Initial window size for opengl display */
#define RANDOM_SEED 1	/* 1 if random initialization, 0 otherwise */
#define SEED 40986	/* 4227 124227 7534298 8522011 525642 93413*/
#define MAX_DISP 1	/* To be defined in the interface */

/* Initial shape of the lattice, not yet totally implemented (see in lattice.c for a piece of code
 */
#define ROUND_LATTICE 0	/* Initial shape of the lattice, not yet totally implemented (see in lattice.c for a piece of code) */

#define INTERFACE_ENERGY 1

/* Variable used in T1 process */
#define T1_PROCESS 1	/* Activates T1 process */
#define CUTOFF_INF 0.05 /* 0.005 */
#define CUTOFF_SUP 0
#define NB_T1_RELAX 5

/* Only used to display the progeny of the center cell (1 to display it)
 */
#define CLONE_DIVISION 1

#define CLONE_DIVISION_ONLY 0 /* Divide only the center cell */


/* Varibles concerning the evolution process
 */
#define SIZEABLE_CELLS 0	/* Possibility of having different sizes of cell */
#define DIVISION_RATE_CHANGEABLE 0	/* Possibility of changing the rate of division */


#define AXIS_DIVISION 1 /* Divide along the main axis of the cell */

/* Pressure killing parameters
 */
#define PRESSURE_KILL 0	/* Activates pressure killing mechanism */
#define PRESSURE_THRESHOLD 0.6	/* Limit pressure (cells whose pressure are higher are killed) */
#define PRESSURE_KILL_INTER 100	/* Interval between 2 killing checks */

/* If true, the displacement due to the bending energy is also saved
 */
#define SAVE_BENDING 0

/* Constants of transverse evolutions
 */
#define THICK_BOUND 1.0

/* Variable used to dislay several quantities (see opengl.c)
 */
#define DISPLAYFORCECOEF 1
#define DISTORTCOEF 1
#define ANTISYMCOEF 1

#define PRESSURECOEF 0.5

#define DIFFER 1

#endif /* CONST_H */
