#ifndef OPTI_H
#define OPTI_H 1

int check_t1();
void update_pa();
//contributions from cell area terms
double area_fct(Dual *pcell);
double area_fctp(Dual *pcell);
double d2E_dA2(Dual *pcell);
int area_exceptions(Dual *pcell);
int area_models(Dual *pcell,double *fct);
//contributions from cell perimeter terms
double rho_fct(Dual *pcell);
double rho_fctpp(Dual *pcell);
double d2E_dP2(Dual *pcell);
int perim_exceptions(Dual *pcell);
int perim_models(Dual *pcell,double *fct);
//contributions from bond terms
double rhol_fct(Bond *pbond);
double rhol_fctpl(Bond *pbond);
double d2E_dL2(Bond *pbond);
int bond_exception(Bond *pbond);
int bond_models(Bond *pbond,double *fct);

//contributions from external potentials
double extern_pot_fct(Dual *pcell);
vector3 extern_pot_fctp(Node *pvert);

//contribution from ``thermal'' fluctuations
double init_thermalfluc(); //finds initial mean energy (in units of force fluctuation) for thermal system
void thermalfluc(gsl_vector *x,double MeanEnergy,double CoeffOfVariation);

//overall energy and denergy computations
double energy(const gsl_vector *x,void *params);
void denergy(const gsl_vector *x, void *params, gsl_vector *df);
void enerdener(const gsl_vector *x, void *params, double *f, gsl_vector *df);
//optimization functions
int optimize(void);
void optimizet1();

// geometric terms when taking derivatives
int dP_dr(Dual *pcell,int idx,double *outvec);
int d2P_drdr(Dual *pcell,int idx1,int idx2,double **outmat);

int parallelProj(Node *pv1,Node *pv2,double *outvec);
int perpProj(Node *pv1,Node *pv2,double **outmat);
int d2RotatePerpProj(Node *pv1,Node *pv2,double *outvec);

int dA_dr(Dual *pcell,int idx,double *outvec);
int d2A_drdr(Dual *pcell,int idx1,int idx2,double **outmat);

int checkModelTypes();

#endif /* OPTI_H */
