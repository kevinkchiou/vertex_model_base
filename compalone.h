#ifndef SRC_COMPALONE_H
#define SRC_COMPALONE_H 1

/* In this file are defined all the preprocessor macros used in
 *  compalone.c. These are constants or switches to activate a process
 */

/* Types of action */
#define CLONE_STAT 0  /* Statistics on clones : nb of cells in clones, pressure, etc */
#define CLONE_STAT_GAMMA 1 /* Statistics when changing the division rate of a clone */
#define CLONE_STAT_GAMMA_DIFF 0
#define CLONE_STAT_KILL g0 /* Statistics with mechanical feedback */
#define CLONE_STAT_SIZE 0 /* Statistics when changing the cell size of a clone */
#define RHEO 0 /* Rheological analysis */
#define SIMULATION_SEQUENCE 0 /* Sequence of simulations */
#define NEWH_STAT 0

/* Preprocessor macros used in CLONE_STAT and CLONE_STAT_GAMMA*/
#define NB_WEB_TOT 100	        /* # of webs to generate */
#define NB_CELL_LIM 300 	/* # of cells in the starting lattice */
#define NB_RELAX 300 		/* nb of steps to relax a division */
#define NB_CELLS_IN_CLONE 1 /* nb of cells in the clone at the beginning */

/* Preprocessor macros used in CLONE_STAT */
#define NB_CLONES_PER_WEB 5 	/* # of clones studied by web */
#define AV_CLONE_LIM 33 	/* Maximum number of cells by clone on average */

/* Preprocessor macros used in CLONE_STAT_GAMMA */
#define GAMMA_CLONE_0 6000	/* Starting value of gamma */
#define GAMMA_CLONE_LIM 6000.1	/* Final value of gamma */
#define GAMMA_CLONE_INC 0.4	/* increment of gamma */
#define NB_CELL_LIM_2 1600	/* Final number of cells */
#define RADIUS_2 4

/* Preprocessor macros used in CLONE_STAT_SIZE */
#define SIZE_CLONE_0 2	/* Starting value of gamma */
#define SIZE_CLONE_LIM 2.1	/* Final value of gamma */
#define SIZE_CLONE_INC 0.4	/* increment of gamma */


/* Preprocessor macros used in/for EVOLUTION_ONLY */
#define EVOLUTION_ONLY 0
#define NB_CELL_LIM_EVOLUTION 2816

/* Data for series of simulations
 */
#define NB_SIMUL_BACK 1
#define NBS_BACK {150}
#define NB_SIMUL_CLONE 1
#define NBS_CLONE {1}
#define NBS_CLONE2 {15}


#endif /* SRC_COMPALONE_H */

/* Linear feedback
 */
#define PRESSURE_MIN  0.0 //0.2
#define PRESSURE_MAX  0.2 //0.7
#define PRESSURE_0 0.0

/* Square feedback
 */
#define PRESSURE_MID  0.0
#define PRESSURE_WTH  0.2

/* Apoptosis
 */
#define PRES_KILL_MIN  -0.1 //0.0
#define PRES_KILL_MAX   0.5 //1.0


extern gsl_rng *rng;

double feedback_coeff(Dual *pcell);
Dual *linear_feedback();

Dual *wo_feedback();
//Dual *wh_feedback();
Dual *pre_kill_cell();
int kill_divide_cell();
double linear_feedback_coeff();
double pres_coeff_kill();
double pres_coeff();
void write_pressure(FILE *ptfile, int tsteps);
void update_areas();
int kill_divide_cell_p(FILE *lofile);
double pressureavg();
int divc_selective();
int divc_stocha();
int divc_feedback();
void mod_tension_relax(double eps);
