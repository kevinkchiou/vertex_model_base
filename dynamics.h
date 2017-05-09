
typedef struct
{
	double delta;
	double notch;
	double reg;
}Lat;

typedef struct
{
    double bn;
    double bd;
    double t;
    double kc;
    double kr;
    double br;
}Parameter;

extern Lat **web_lat;

Parameter *allocparam();

Lat **createlatvector(int dim);
Lat **realloclatvector(Lat **latvector,int dim);
Lat *alloclat();

int averageND(const Dual *pcell,const double y[],double *avgN,double *avgD);
int func(double t, const double y[], double f[], void *params);
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params);
void dynamics();
void init_dynamics();
void init_dynamics_rnd(int initialize);
void setclone_dynamics();
void division_dynamics(int origcellnum);
int hstate_check();
