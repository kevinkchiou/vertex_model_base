
void coord_mod();//modifies coords by an offset in z

void border_check();//checks to see which vertices are on the border for opti
void division_cell_slow(Dual *pcell, int ivertex);//other division algorithm to maintain planarity

void tension_modify();
void tension_modify1();
double tension_modify2(double eps);
void tension_uniform();
void tension_uniform_eraseclone();
void kagome_pre_modify();
void freeparking(int **park);
int **rnd_parking(int *oc, int *unoc, int *oclen,int *unoclen);
void kagome_nn_modify();
void kagome_nnn_modify();
void kagome_modify2();
void kagome_dblmodify();
void kagome_reg();
void init_lattice2(int *array, int num);
void import_lattice(char * filename);
void kagomet1_dd();
void kagomet1_nd();
void pulse_event(double step_size,double noise_strength);
int *division_select_once(int list_cells_to_divide[],int divlistlen,int *pnum);
void torus_division_set_parameters();

Dual *findCenterCell();

//following functions are used to create the dynamical matrix
gsl_matrix* computeDynMat();
//outerProd() is used for quadratic terms of Hamiltonian (ie.
//the outer products of sigma_ij = dP_i dP_j, dA_i * dA_j
int outerProd(double *vec1,double *vec2,double **outmat);
//addIntoDynMat() adds individual terms from smaller matrix
//inmat into the appropriate locations in dynmat
int addIntoDynMat(gsl_matrix *dynmat,double **inmat,int idx1,int idx2,double prefactor);

int issameMatrices(const gsl_matrix *mat1,const gsl_matrix *mat2,int dim1,int dim2);

//following functions create compatibility matrices
int addVal_gsl_matrix(gsl_matrix *m,size_t i,size_t j,double val);
gsl_matrix* computePerimCompatMat();
