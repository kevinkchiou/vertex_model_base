
typedef struct
{
	double rho;
}Evolparam;

Evolparam *allocevolparam();

//functions for performing gradient descent evolution
//requires model-dependent functions from opti.c
int evolfunc(double t, const double y[], double f[], void *params);
int evoljac(double t, const double y[], double *dfdy, double dfdt[], void *params);
int hessian_length2d(gsl_matrix *m,int vertexindex,const double y[]);
int hessian_length3d(gsl_matrix *m,int vertexindex,const double y[]);
int hessian_area2d(gsl_matrix *m,int cellindex,const double y[]);
int hessian_area3d(gsl_matrix *m,int cellindex,const double y[]);
int test_symmetric(gsl_matrix *m);
int test_eig(gsl_matrix *m,gsl_vector *eval,gsl_matrix *evec);

void evol(const double delta_t,const char *method);

//functions for stress-activated contraction
void evol_stress_activation(const double delta_t);
void lead_contract(int *celllist,int celllist_length);
void init_evol_stress_activation();
void check_cell_state_stress_activation(double baseline,double threshold,double dt);

//functions for simultaneous evolution
void init_evol_total_activation(Dual *pcell);
void evol_total_activation(const double dt,const double E0);

//introducing positional fluctuations
void posfluc(double CoeffofVariation);
