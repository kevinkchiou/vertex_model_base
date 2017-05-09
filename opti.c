#include <string.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "definelat.h"
#include "lattice.h"
#include "locerror.h"
#include "save.h"
#include "measurements.h"
#include "pngwrite.h"
#include "compalone.h"
#include "opti.h"
#include "in_development.h"

double prefactor1 = 100.0;
double prefactor2 = 0.0;
double radius = 8.0;
/* Make a t1 process if necessary
 */
int check_t1()
{
   int i,j;
   static int nbt1=0;
   Node *pvert, *pvertn;

   for (i=0;i<nb_vertex_tot;i++){
      pvert=web[i];
      for (j=0;j<3;j++){
         pvertn=pvert->pneighb[j];
         if (dist(pvert,pvertn)<CUTOFF_INF && flagt1<4){
            printf("performing t1 transform on vertex %d and neighbor %d\n",indexvertexinweb(pvert),j);
            t1flagcheck(indexvertexinweb(pvert),j);//increments flagt1 when necessary to prevent loop
            t1transform(pvert,j);
            //border_check();
            nbt1++;
            return 1;
         }
      }
   }
   if (flagt1>=4){flagt1=0;}
   return 0;
}

/* update perimeter and area
 */
void update_pa()
{
    int i;
    double carea;
    Dual *pcell;
    for (i=0;i<nb_cell_tot;i++){
	pcell=web_dual[i];
	carea=area(pcell);
	pcell->area=carea;
	pcell->sqrtarea=sqrt(carea);
	pcell->perimeter=perimeter(pcell);
	pcell->centroid=centroid(pcell);
    }
}

// area term contributions to energy and its derivatives
double area_fct(Dual *pcell){
    double fct[3];
	//consistent way of computing area models contributions. Get zeroth derivative
	if(area_models(pcell,fct)==0){return 0.;}
    return fct[0];
}

double area_fctp(Dual *pcell){
    double fct[3];
	//consistent way of computing area models. Get first derivative
	if(area_models(pcell,fct)==0){return 0.;}
    return fct[1];
}

double d2E_dA2(Dual *pcell){
	double fct[3];
	//Consistent way of computing area models. Get second derivative
	if(area_models(pcell,fct)==0){return 0.;}
	return fct[2];
}

int area_exceptions(Dual *pcell){
	//write exceptions to area contributions
	if(pcell->marker.border==1){return 1;}
	if(TOROIDAL==0 && pcell->idx==0){return 1;}

	return 0;
}

int area_models(Dual *pcell,double *fct)
{
	int type;
	double A,Asoll,kA;
	if(strcmp(AREAMODELTYPE,"none")==0 || area_exceptions(pcell)==1 || fabs(pcell->kA)<1e-9){fct[0]=0.;fct[1]=0.;fct[2]=0.;return 0;}
	else if(strcmp(AREAMODELTYPE,"quadratic")==0){type=1;A=pcell->area;Asoll=pcell->area_soll,kA=pcell->kA;}
	else if(strcmp(AREAMODELTYPE,"hyperelastic")==0){type=2;A=pcell->area;Asoll=pcell->area_soll;kA=pcell->kA;}
	else locerror("area_models()","Not a valid area model type!");

	switch(type) {
		case 1 :
			fct[0]=0.5*kA*(A-Asoll)*(A-Asoll);
			fct[1]=kA*(A-Asoll);
			fct[2]=kA;		
			return 1;
		case 2 :
			fct[0]=0.5*kA*(A-Asoll)*(A-Asoll)/sqrt(A);
			fct[1]=kA*((A-Asoll)/sqrt(A) - 0.25*(A-Asoll)*(A-Asoll)/sqrt(A*A*A));
			fct[2]=(kA/sqrt(A))*(Asoll/A+3/8*(1.-Asoll/A)*(A-Asoll/A));
			return 2;
	}
}

// perimeter term contributions to energy and derivatives
double rho_fct(Dual *pcell){
	double fct[3];
	if(perim_models(pcell,fct)==0){return 0.;}
	return fct[0];
}

double rho_fctpp(Dual *pcell){
	double fct[3];
	if(perim_models(pcell,fct)==0){return 0.;}
	return fct[1];
}

double d2E_dP2(Dual *pcell){
	double fct[3];
	if(perim_models(pcell,fct)==0){return 0.;}
	return fct[2];
}

int perim_exception(Dual *pcell){
	if(pcell->marker.border==1){return 1;}
	return 0;
}

int perim_models(Dual *pcell,double *fct){
	int type;
	double P,Psoll,kP,sigma;
	if(perim_exception(pcell)==1 || strcmp(PERIMMODELTYPE,"none")==0){fct[0]=0.;fct[1]=0.;fct[2]=0.;return 0;}
	else if(strcmp(PERIMMODELTYPE,"linear")==0){type=1;P=pcell->perimeter;sigma=pcell->sigma;}
	else if(strcmp(PERIMMODELTYPE,"quadratic")==0){type=2;P=pcell->perimeter;Psoll=pcell->perim_soll;kP=pcell->kP;}

	switch(type) {
		case 1 :
			fct[0]=sigma*P;
			fct[1]=sigma;
			fct[2]=0.;
			return 1;
		case 2 :
			fct[0]=0.5*kP*(P-Psoll)*(P-Psoll)+sigma*P;
			fct[1]=kP*(P-Psoll)+sigma;
			fct[2]=kP;
			return 2;
	}
}

// Bond contributions to energy and derivatives
double rhol_fct(Bond *pbond)
{
	double fct[3];
	if(bond_models(pbond,fct)==0){return 0.;}
	return fct[0];
}

double rhol_fctpl(Bond *pbond)
{
	double fct[3];
	if(bond_models(pbond,fct)==0){return 0.;}
	return fct[1];
}

double d2E_dL2(Bond *pbond)
{
	double fct[3];
	if(bond_models(pbond,fct)==0){return 0.;}
	return fct[2];
}

int bond_exceptions(Bond *pbond){
	if(fabs(pbond->kappa)<1e-9 && fabs(pbond->lambda)<1e-9){return 1;}
	return 0; //no exceptions
}

int bond_models(Bond *pbond,double *fct){
	double kappa,L0,lambda,L;
	if(bond_exceptions(pbond)==1){fct[0]=0.;fct[1]=0.;fct[2]=0.;return 0;}
	else{kappa=pbond->kappa;lambda=pbond->lambda;L0=pbond->L0;L=dist(pbond->pnvert[0],pbond->pnvert[1]);}
	
	//so far, just the one model
	fct[0]=0.5*kappa*(L-L0)*(L-L0)-lambda*L;
	fct[1]=kappa*(L-L0)-lambda;
	fct[2]=kappa;
	return 1;
}

// contributions from external potentials
double extern_pot_fct(Dual *pcell){
    double ext=0.0,x,y,z;
    int i,num = pcell->nb_vertices;
    Node *pvert;

    ext = 0.0;
    /*
       for(i=0;i<num;i++){
       pvert = pcell->vertexlist[i];
       x = pvert->x;y = pvert->y;z = pvert->z;
       ext+=0.5*prefactor1*(sqrt(x*x+y*y+z*z)-radius)*(sqrt(x*x+y*y+z*z)-radius);
       }*/
    return ext;//(distance-radius)^2 is the parameter
}

vector3 extern_pot_fctp(Node *pvert){
    vector3 dpot;
    double x,y,z;

    dpot.x=dpot.y=dpot.z=0.0;//for flat external potential case
    /*
       x=pvert->x;y=pvert->y;z=pvert->z;

       dpot.x=prefactor1*(sqrt(x*x+y*y+z*z)-radius);
       dpot.y=prefactor1*(sqrt(x*x+y*y+z*z)-radius);
       dpot.z=prefactor1*(sqrt(x*x+y*y+z*z)-radius);
     */
    return dpot;
}

/* ------------------------------------------------------------------------------ */
//``thermal'' fluctuations
double init_thermalfluc(){
	int i;
	double E0,F0,L0=0.;
	Dual *pcell=NULL;
	gsl_vector *x;
	void *params;

	x=gsl_vector_alloc(3*nb_vertex_tot);
	for(i=0;i<nb_vertex_tot;i++){
		gsl_vector_set(x,3*i,web[i]->x);
		gsl_vector_set(x,3*i+1,web[i]->y);
		gsl_vector_set(x,3*i+2,web[i]->z);
	}

	E0=energy(x,params); //total energy in the system
	//find the appropriate forces by finding relevant lengthscales
	for(i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		L0+=(pcell->perimeter)/2.; //total ``length'' in the system
	}
	
	F0=E0/L0*((double)nb_cell_tot/((double)nb_vertex_tot)); //normalized ``mean force'' on each vertex
	gsl_vector_free(x);
	return F0; //normalized ``mean force'' in the system
}

void thermalfluc(gsl_vector *df,double F0,double coeff){
	double deltaF,sigma;
	double temp=0.;
	int i=0;
	
	sigma=coeff*F0;
	for(i=0;i<NUMDIM*nb_vertex_tot;i++){
		deltaF = gsl_ran_gaussian(rng,sigma);
		temp=gsl_vector_get(df,i);
		gsl_vector_set(df,i,temp+deltaF);
	}
}


/* ------------------------------------------------------------------------------ */

double energy(const gsl_vector *x,void *params){   /* Compute the energy */
    int i;
    double h=0;
    Dual *pcell;
	Bond *pbond;

    for (i=0;i<nb_vertex_tot;i++){
	web[i]->x=gsl_vector_get(x,3*i);
	web[i]->y=gsl_vector_get(x,3*i+1);
	web[i]->z=gsl_vector_get(x,3*i+2);
    }
    update_pa();
	
    for (i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		h+=rho_fct(pcell);
		h+=area_fct(pcell)+extern_pot_fct(pcell);
    }
	
	for(i=0;i<nb_bond_tot;i++){
		pbond=web_bond[i];
		h+=rhol_fct(pbond);
	}

    return h;
}

void denergy(const gsl_vector *x, void *params, gsl_vector *df){
    /* Compute the gradient of the energy */
    int i,j;
    double vertx,verty,vertz;
    double dx, dy, dz, dxa, dya, dza, dxd, dyd, dzd;
    double locdist, prefact;
    double ax,ay,az,bx,by,bz;
    double xcross,ycross,zcross;
    Node *pvert, *pvert_prev, *pvert_next;
    Dual *pcell0, *pcell;
	Bond *pbond_prev,*pbond_next;
    vector3 cent,dpot;

    for (i=0;i<nb_vertex_tot;i++){
	web[i]->x=gsl_vector_get(x,3*i);
	web[i]->y=gsl_vector_get(x,3*i+1);
	web[i]->z=gsl_vector_get(x,3*i+2);
    }
    update_pa();
    pcell0=web_dual[0];
    for (i=0;i<nb_vertex_tot;i++){
	pvert=web[i];
	dx=dy=dz=dxa=dya=dza=dxd=dyd=dzd=0;
	for (j=0; j<3; j++){
	    pcell    = pvert -> pncell[j];
		pbond_next = pvert->pnbond[j];
	    pvert_next = pvert ->pneighb[j];
		pbond_prev = pvert->pnbond[(j+1)%3];
	    pvert_prev = pvert -> pneighb[(j+1)%3];
	    cent = pcell->centroid;

		//check of stuff!
		//if(pbond_next->pnvert[0]!=pvert_next && pbond_next->pnvert[1]!=pvert_next){locerror("opti","bond error");}
		//if(pbond_prev->pnvert[0]!=pvert_prev && pbond_next->pnvert[1]!=pvert_prev){locerror("opti","bond error");}

	    // Displacement along x, y, z due to perimeter extension 
		locdist=dist(pvert_next,pvert);
		prefact = rho_fctpp(pcell)/locdist;
		prefact += 0.5*rhol_fctpl(pbond_next)/locdist;
		dx += distx(pvert,pvert_next)*prefact;
		dy += disty(pvert,pvert_next)*prefact;
		dz += distz(pvert,pvert_next)*prefact;

		locdist=dist(pvert_prev,pvert);
		prefact = rho_fctpp(pcell)/locdist;
		prefact += 0.5*rhol_fctpl(pbond_prev)/locdist;
		dx += distx(pvert,pvert_prev)*prefact;
		dy += disty(pvert,pvert_prev)*prefact;
		dz += distz(pvert,pvert_prev)*prefact;


		// Displacement along x and y due to area extension 

		ax=distx(pvert_next,pvert);ay=disty(pvert_next,pvert);az=distz(pvert_next,pvert); 
		bx=centdistx(cent,pvert);by=centdisty(cent,pvert);bz=centdistz(cent,pvert);
		zcross = ax*by-ay*bx; xcross = ay*bz-az*by; ycross = az*bx-ax*bz;
		prefact = area_fctp(pcell); //CHECK_HERE as well
		//dxa += prefact*0.5*(pvert_next -> y - pvert_prev -> y); 
		//dya += prefact*0.5*(pvert_prev -> x - pvert_next -> x);
		//dza += prefact*0.5*(pvert_prev -> z - pvert_next -> z);
		dxa+=prefact*0.5*(zcross*(ay-by)-ycross*(az-bz))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dya+=prefact*0.5*(xcross*(az-bz)-zcross*(ax-bx))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dza+=prefact*0.5*(ycross*(ax-bx)-xcross*(ay-by))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);

		ax=centdistx(cent,pvert);ay=centdisty(cent,pvert);az=centdistz(cent,pvert); 
		bx=distx(pvert_prev,pvert);by=disty(pvert_prev,pvert);bz=distz(pvert_prev,pvert);
		zcross = ax*by-ay*bx; xcross = ay*bz-az*by; ycross = az*bx-ax*bz;
		dxa+=prefact*0.5*(zcross*(ay-by)-ycross*(az-bz))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dya+=prefact*0.5*(xcross*(az-bz)-zcross*(ax-bx))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dza+=prefact*0.5*(ycross*(ax-bx)-xcross*(ay-by))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		

	}

	dpot=extern_pot_fctp(pvert);
	dxd+=dpot.x;
	dyd+=dpot.y;
	dzd+=dpot.z;

	gsl_vector_set(df, 3*i, dx+dxa+dxd);
	gsl_vector_set(df, 3*i+1, dy+dya+dyd);
	gsl_vector_set(df, 3*i+2, dz+dza+dzd);

	forcevector[0][i]=gsl_vector_get(df,3*i);
	forcevector[1][i]=gsl_vector_get(df,3*i+1);
	forcevector[2][i]=gsl_vector_get(df,3*i+2);
	}
}


void enerdener(const gsl_vector *x,void *params,double *f, gsl_vector *df){

	*f=energy(x,params);
	denergy(x,params,df);
}

int optimize(void) {
    int i;
    int iter=0;
    int status;
    gsl_vector *startpos;
    gsl_multimin_function_fdf hamilton;
    const gsl_multimin_fdfminimizer_type *minitype;
    gsl_multimin_fdfminimizer *minimize;
    double par[1] = {(double) MUP};

    int N=3*nb_vertex_tot;
    hamilton.f      = &energy;
    hamilton.df     = &denergy;
    hamilton.fdf    = &enerdener;
    hamilton.n      = N;
    hamilton.params = &par;


	printf("nb_bond_tot = %d\n",nb_bond_tot);fflush(stdout);
    startpos = gsl_vector_alloc(N);
    for (i=0;i<nb_vertex_tot;i++){
	gsl_vector_set(startpos,3*(i),web[i]->x);
	gsl_vector_set(startpos,3*(i)+1,web[i]->y);
	gsl_vector_set(startpos,3*(i)+2,web[i]->z);
    }
    minitype=gsl_multimin_fdfminimizer_conjugate_pr;
    //   minitype=gsl_multimin_fdfminimizer_vector_bfgs;
    //   minitype=gsl_multimin_fdfminimizer_conjugate_fr;
    minimize=gsl_multimin_fdfminimizer_alloc(minitype, N);

    gsl_multimin_fdfminimizer_set(minimize, &hamilton, startpos, 0.00001, 0.00001);
    do {
	iter++;
	status = gsl_multimin_fdfminimizer_iterate(minimize);

   //if(iter==500){printf("iter=500\n");}
   //if(iter==750){printf("iter=750\n");}
   if(iter==999){printf("iter=999, it is about to break...\n");}
	if (status) break;

	status = gsl_multimin_test_gradient(minimize->gradient, 1e-4);
    } while (status==GSL_CONTINUE && iter<1000);
    for (i=0;i<nb_vertex_tot;i++){
	web[i]->x=gsl_vector_get(minimize->x,3*(i));
	web[i]->y=gsl_vector_get(minimize->x,3*(i)+1);
	web[i]->z=gsl_vector_get(minimize->x,3*(i)+2);
    }

    gsl_multimin_fdfminimizer_free(minimize);
    gsl_vector_free(startpos);
    return 0;
}

void optimizet1(){
    int icheck=1; 
    optimize();
    do {
	if (check_t1()){
	    optimize();
	    icheck=1;
	}
	else icheck=0;
    } while (icheck);
}

int dP_dr(Dual *pcell,int idx,double *outvec){
	Node *pv,*pvn,*pvp;
	double *temp;
	int i;

	temp=(double*)calloc(NUMDIM,sizeof(double));

	pv=pcell->vertexlist[idx];
	pvn=pcell->vertexlist[(idx+1)%(pcell->nb_vertices)];
	pvp=pcell->vertexlist[(idx-1)%(pcell->nb_vertices)];

	parallelProj(pv,pvn,temp);
	for(i=0;i<NUMDIM;i++){outvec[i]=temp[i];}
	parallelProj(pv,pvp,temp);
	for(i=0;i<NUMDIM;i++){outvec[i]+=temp[i];}
	
	free(temp);
	return 0;
}

int d2P_drdr(Dual *pcell,int idx1,int idx2,double **outmat){
	int i,j,flag;
	Node *pv1,*pv2;
	double **temp;
	int idx;

	//Pre-stress prefactor for this second order term
	//prefac=rho_fctpp(pcell);
	//// NO LONGER DO THIS! This term is meant to be purely geometric now!

	//test to see if the vertices share a bond within this paricular cell
	//if so, compute the contribution and return 1. If not, perform no
	//computations and return 0.
	if(idx1==idx2){
		//allocate memory for computations
		temp=(double**)calloc(NUMDIM,sizeof(double*));
		for(i=0;i<NUMDIM;i++){temp[i]=(double*)calloc(NUMDIM,sizeof(double));}

		pv1=pcell->vertexlist[idx];
		pv2=pcell->vertexlist[(idx+1)%pcell->nb_vertices];
		perpProj(pv1,pv2,temp);
		for(i=0;i<NUMDIM;i++){for(j=0;j<NUMDIM;j++){outmat[i][j]=temp[i][j];}}
		pv2=pcell->vertexlist[(idx-1)%pcell->nb_vertices];
		perpProj(pv1,pv2,temp);
		for(i=0;i<NUMDIM;i++){for(j=0;j<NUMDIM;j++){outmat[i][j]+=temp[i][j];}}
		flag=1;

		free(temp);//free memory used for computation
	}
	else if(idx1==((idx2+1)%pcell->nb_vertices) || idx1==((idx2-1)%pcell->nb_vertices)){
		//allocate memory for computations
		temp=(double**)calloc(NUMDIM,sizeof(double*));
		for(i=0;i<NUMDIM;i++){temp[i]=(double*)calloc(NUMDIM,sizeof(double));}

		pv1=pcell->vertexlist[idx1];pv2=pcell->vertexlist[idx2];
		perpProj(pv1,pv2,temp);
		for(i=0;i<NUMDIM;i++){for(j=0;j<NUMDIM;j++){outmat[i][j]=-temp[i][j];}}
		flag=1;

		free(temp);//free memory used for computation
	}
	else{flag=0;}
	return flag;
}

int parallelProj(Node *pv1,Node *pv2,double *outvec){
	int i,idx=-1;
	double dx,dy,dz,d;
	for(i=0;i<3;i++){if(pv1->pneighb[i]==pv2){idx=i;}}
	if(idx==-1){locerror("parallelProj()","Designated vertices do not share a bond!\n");return 0;}
	dx=distx(pv1,pv2);dy=disty(pv1,pv2);
	if(NUMDIM==3){dz=distz(pv1,pv2);}
	else if(NUMDIM==2){dz=0.;}
	d=sqrt(dx*dx+dy*dy+dz*dz);

	outvec[0]=dx/d; //x-component
	outvec[1]=dy/d; //y-component
	if(NUMDIM==3) outvec[2]=dz/d; //z-component

	return 1;
}

int d2RotatePerpProj(Node *pv1,Node *pv2,double *outvec){
	int i,idx=-1;
	double dx,dy,d;
	for(i=0;i<3;i++){if(pv1->pneighb[i]==pv2){idx=i;}}
	if(idx==-1){locerror("d2RotatePerpProj","Designated vertices do not share a bond!\n");return 0;}
	dx=-1.*disty(pv1,pv2);dy=distx(pv1,pv2); //only in two dimensions, pi/2 rotation ccw
	d=sqrt(dx*dx+dy*dy);

	outvec[0]=dx/d;outvec[1]=dy/d;
	return 1;
}

int perpProj(Node *pv1,Node *pv2,double **outmat){
	int i,idx=-1;
	double dx,dy,dz,d2,d;

	for(i=0;i<3;i++){if(pv1->pneighb[i]==pv2){idx=i;}}
	if(idx==-1){locerror("perpProj()","Designated vertices do not share a bond!\n");}
	dx=distx(pv1,pv2);dy=disty(pv1,pv2);
	if(NUMDIM==3){dz=distz(pv1,pv2);}
	else if(NUMDIM==2){dz=0.;}
	d2=dx*dx+dy*dy+dz*dz;d=sqrt(d2);

	outmat[0][0]=(1.-dx*dx/d2)/d;
	outmat[0][1]=(-dx*dy/d2)/d;
	outmat[1][0]=outmat[0][1];
	outmat[1][1]=(1.-dy*dy/d2)/d;
	if(NUMDIM==3){
		outmat[0][2]=(-dx*dz/d2)/d;
		outmat[2][0]=outmat[0][2];
		outmat[1][2]=(-dy*dz/d2)/d;
		outmat[2][1]=outmat[1][2];
		outmat[2][2]=(1.-dz*dz/d2)/d;
	}
	return 1;
}

int dA_dr(Dual *pcell,int idx,double *outvec){
	Node *pvn,*pvp;
	int i;

	pvn=pcell->vertexlist[(idx+1)%(pcell->nb_vertices)];
	pvp=pcell->vertexlist[(idx-1)%(pcell->nb_vertices)];

	if(NUMDIM==2){
		//this comes from r_p cross r, and r cross r_n
		//giving positive areas (ie ccw orientation)
		outvec[0]=0.5*disty(pvp,pvn);
		outvec[1]=-0.5*distx(pvp,pvn);
	}
	else if(NUMDIM==3){ //much more complicated, reduces to the d=2 case
		//should match area extension terms in the denergy() function
		//main difference is that A = \sqrt{\vec{A} \cdot \vec{A}} where
		//\vec{A} results from the d=3 cross product of node positions
	}
	return 1;
}

int d2A_drdr(Dual *pcell,int idx1,int idx2,double **outmat){
	Node *pv1,*pv2;
	double fac;
	
	if(idx1!=((idx2+1)%pcell->nb_vertices) && idx1!=((idx2-1)%pcell->nb_vertices)){
		//if they are not neighboring indices in this particular cell, return 0
		//otherwise the function will continue to compute the terms
		return 0;
	}
	//For d=2, no geometric information is needed for this term, but is 
	//important for d=3
	pv1=pcell->vertexlist[idx1];pv2=pcell->vertexlist[idx2];
	fac=(pcell->area - pcell->area_soll);
	if(NUMDIM==2){
		outmat[0][0]=0.;outmat[1][1]=0.;
		if(idx1<idx2){ //check orientation
			outmat[0][1]=0.5*fac;
			outmat[1][0]=-0.5*fac;
		}
		else if(idx1>idx2){
			outmat[0][1]=-0.5*fac;
			outmat[1][0]=0.5*fac;
		}
		else{locerror("d2A_drdr","Indices cannot be the same!\n");}

	}
	else if(NUMDIM==3){ //tremendously more complicated! 
		//in d=3 there are terms like 1/|A| \vec{A} \cdot 
		// \frac{\partial{\vec{A}}}{\partial{\vec{r}}}. 
		//roughly all the terms here include
		//1/(2*A^3)*(A\cdot\partial_i{A})*(A\cdot\partial_j{A}) +
		//1/A*(\partial_i{A}\cdot\partial_j{A} + A\cdot\partial_{ij}{A})
	}
	return 1;
}

int checkModelTypes(){
	char a[50],p[50];
	strncpy(a,AREAMODELTYPE,50);
	strncpy(p,PERIMMODELTYPE,50);
	if(strcmp(a,"quadratic")!=0 && strcmp(a,"none")!=0 && strcmp(a,"hyperelastic")!=0){
		printf("Area energy model type = %s\n",a);
		locerror("checkModelTypes()","Not a valid area energy model!");
	}
	if(strcmp(p,"quadratic")!=0 && strcmp(p,"linear")!=0 && strcmp(p,"none")!=0){
		printf("Perimeter energy model type = %s\n",p);
		locerror("checkModelTypes()","Not a valid perimeter energy model!");
	}
	return 0;
}
