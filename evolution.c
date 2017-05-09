#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>

#include "definelat.h"
#include "locerror.h"
#include "lattice.h"
#include "measurements.h"
#include "evolution.h"
#include "compalone.h"
#include "opti.h"



Evolparam *allocevolparam(){
    Evolparam *a;
    a=(Evolparam *)malloc(sizeof(Evolparam));
    return a;
}

int evolfunc(double t, const double y[], double f[], void *par){
    int i,j,k,nbneighb;
    Node *pvert;
    double *result;
	gsl_vector *x,*df;
	int numdim = NUMDIM;
	double E0;

	E0=((double*)par)[0];

	x  = gsl_vector_alloc(3*nb_vertex_tot);
	df = gsl_vector_alloc(3*nb_vertex_tot);
	//computes elastic (hamiltonian)  contributions
    for(i=0;i<nb_vertex_tot;i++){
		gsl_vector_set(x,3*i,y[numdim*i]);
		gsl_vector_set(x,3*i+1,y[numdim*i+1]);
		if(numdim==3){gsl_vector_set(x,3*i+2,y[numdim*i+2]);}
		else{gsl_vector_set(x,3*i+2,0.0);}
	}
	denergy(x,par,df); //use opti.c component
	//thermalfluc(df,E0,0.2); //add in thermal fluctuations to df vector
	for(i=0;i<nb_vertex_tot;i++){
		f[numdim*i]=-1.0*gsl_vector_get(df,3*i);
		f[numdim*i+1]=-1.0*gsl_vector_get(df,3*i+1);
		if(numdim==3){f[numdim*i+2]=-1.0*gsl_vector_get(df,3*i+2);}
   		forcevector[0][i]=gsl_vector_get(df,3*i);
   		forcevector[1][i]=gsl_vector_get(df,3*i+1);
   		forcevector[2][i]=gsl_vector_get(df,3*i+2);
	}

	gsl_vector_free(x);gsl_vector_free(df);
    return GSL_SUCCESS;
}

int evoljac(double t, const double y[], double *dfdy, double dfdt[], void *par){
    int i,j,k,nbneighb,ncell1,ncell2;
    Node *pvert;
    double *result;
	int numdim=NUMDIM;
	int test=0;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,numdim*nb_vertex_tot,numdim*nb_vertex_tot);
    gsl_matrix *m = &dfdy_mat.matrix;
	gsl_matrix_set_zero(m);

    for(i=0;i<nb_vertex_tot;i++){
		if(numdim==2){hessian_length2d(m,i,y);}
		if(numdim==3){hessian_length3d(m,i,y);}
        dfdt[numdim*i]=0.0;dfdt[numdim*i+1]=0.0;
		if(numdim==3){dfdt[numdim*i+2]=0.0;}
	}
	for(i=0;i<nb_cell_tot;i++){
		if(numdim==2){hessian_area2d(m,i,y);}
		if(numdim==3){hessian_area3d(m,i,y);}
    }

	//for testing purposes
	if(test==1){
		j=0;k=0;
		gsl_vector *eval = gsl_vector_alloc(m->size1);
		gsl_matrix *evec = gsl_matrix_alloc(m->size1,m->size2);
		test_eig(m,eval,evec);
		//stuff to test

		for(i=0;i<eval->size;i++){
			if(gsl_vector_get(eval,i)>0.0001){j++;
				printf("Positive eval = %lf...\n",gsl_vector_get(eval,i));fflush(stdout);
			}
			if(gsl_vector_get(eval,i)<-0.0001){k++;
				//printf("Negative eval = %lf...\n",gsl_vector_get(eval,i));fflush(stdout);
			}
		}
		//printf("There are %d positive modes and %d negative modes!\n",j,k);fflush(stdout);

		gsl_vector_free(eval);gsl_matrix_free(evec);

	}

    return GSL_SUCCESS;
}

int hessian_length2d(gsl_matrix *m,int vertindex,const double y[]){

	Node *pvert,*pnvert;
	int i,j,nbneighb,nidx;
	double sig;
	double dx,dy,d,d2;
	double dxdx_d2,dxdy_d2,dydy_d2;
	
	i=vertindex;
	pvert=web[i];

	nbneighb=3;
	for(j=0;j<nbneighb;j++){
		pnvert = pvert->pneighb[j];
		nidx = pnvert->idx;
		dx = ptptdistx(y[2*i],y[2*nidx]);
		dy = ptptdisty(y[2*i+1],y[2*nidx+1]);
		d = sqrt(dx*dx+dy*dy);
		d2=d*d;



		/////// dH_dL * d2L_dr2 terms ////////
		sig = -1.0*tension(pvert->pnbond[j])/d; //tension = dH_dL
		dxdx_d2=dx*dx/d2;dxdy_d2=dx*dy/d2;dydy_d2=dy*dy/d2;

		//i==j terms (need sums from each neighboring vertex)
		m->data[2*i*m->tda+(2*i)]+=sig*(1-dxdx_d2);
		m->data[2*i*m->tda+(2*i+1)]-=sig*dxdy_d2;
		m->data[(2*i+1)*m->tda+2*i]-=sig*dxdy_d2;
		m->data[(2*i+1)*m->tda+(2*i+1)]+=sig*(1-dydy_d2);

		//i!=j terms
		m->data[2*i*m->tda+(2*j)]-=sig*(1-dxdx_d2);
		m->data[2*i*m->tda+(2*j+1)]+=sig*dxdy_d2;
		m->data[(2*i+1)*m->tda+2*j]+=sig*dxdy_d2;
		m->data[(2*i+1)*m->tda+(2*j+1)]-=sig*(1-dydy_d2);


		//////// d2H_dL2 *(dL_dr)^2 terms /////////
		if(d2E_dL2(pvert->pnbond[j])>0.0001){
			//do stuff if this term is non-zero

		}

	}

	return 1;
}

int hessian_length3d(gsl_matrix *m,int vertindex,const double y[]){

	Node *pvert,*pnvert;
	int i,j,nbneighb,nidx;
	double sig;
	double dx,dy,dz,d,d2;
	double dxdx_d2,dxdy_d2,dxdz_d2,dydy_d2,dydz_d2,dzdz_d2;
	
	i=vertindex;
	pvert=web[i];

	nbneighb=3;
	for(j=0;j<nbneighb;j++){
		pnvert = pvert->pneighb[j];
		nidx = pnvert->idx;
		dx = ptptdistx(y[3*i],y[3*nidx]);
		dy = ptptdisty(y[3*i+1],y[3*nidx+1]);
		dz = ptptdistz(y[3*i+2],y[3*nidx+2]);
		d = sqrt(dx*dx+dy*dy+dz*dz);
		d2=d*d;



		/////// dH_dL * d2L_dr2 terms ////////
		sig = -1.0*tension(pvert->pnbond[j])/d; //tension = dH_dL
		dxdx_d2=dx*dx/d2;dxdy_d2=dx*dy/d2;dxdz_d2=dx*dz/d2;
		dydy_d2=dy*dy/d2;dydz_d2=dy*dz/d2;dzdz_d2=dz*dz/d2;

		//i==j terms (need sums from other vertices)
		m->data[3*i*m->tda+(3*i)]+=sig*(1-dxdx_d2);
		m->data[3*i*m->tda+(3*i+1)]-=sig*dxdy_d2;
		m->data[3*i*m->tda+(3*i+2)]-=sig*dxdz_d2;
		m->data[(3*i+1)*m->tda+3*i]-=sig*dxdy_d2;
		m->data[(3*i+1)*m->tda+(3*i+1)]+=sig*(1-dydy_d2);
		m->data[(3*i+1)*m->tda+(3*i+2)]-=sig*dydz_d2;
		m->data[(3*i+2)*m->tda+3*i]-=sig*dxdz_d2;
		m->data[(3*i+2)*m->tda+(3*i+1)]-=sig*dydz_d2;
		m->data[(3*i+2)*m->tda+(3*i+2)]+=sig*(1-dzdz_d2);

		//i!=j terms
		m->data[3*i*m->tda+(3*j)]-=sig*(1-dxdx_d2);
		m->data[3*i*m->tda+(3*j+1)]+=sig*dxdy_d2;
		m->data[3*i*m->tda+(3*j+2)]+=sig*dxdz_d2;
		m->data[(3*i+1)*m->tda+3*j]+=sig*dxdy_d2;
		m->data[(3*i+1)*m->tda+(3*j+1)]-=sig*(1-dydy_d2);
		m->data[(3*i+1)*m->tda+(3*j+2)]+=sig*dydz_d2;
		m->data[(3*i+2)*m->tda+3*j]+=sig*dxdz_d2;
		m->data[(3*i+2)*m->tda+(3*j+1)]+=sig*dydz_d2;
		m->data[(3*i+2)*m->tda+(3*j+2)]-=sig*(1-dzdz_d2);


		//////// d2H_dL2 *(dL_dr)^2 terms /////////
		if(d2E_dL2(pvert->pnbond[j])>0.0001){
			//do stuff if this term is non-zero

		}

	}

	return 1;
}

int hessian_area2d(gsl_matrix *m,int cellindex,const double y[]){

	int a,nbneighb,i,j;
	int idx,idxn,idxp,idxo,idxop,idxon;
	double dx,dy,dxo,dyo;
	double rho,d2E;
	Dual *pcell;
	Node *pv,*pvprev,*pvnext;
	Node *pvo,*pvoprev,*pvonext;

	a = cellindex;
	pcell=web_dual[a];
	rho = pressure(pcell);
	d2E = -1.0*d2E_dA2(pcell);
	nbneighb = pcell->nb_vertices;
	//printf("cell = %d initial declarations finished!\n",a);fflush(stdout);

	for(i=0;i<nbneighb;i++){
		pv=pcell->vertexlist[i];
		pvprev=pcell->vertexlist[((i+nbneighb-1)%nbneighb)];
		pvnext=pcell->vertexlist[((i+1)%nbneighb)];
		idx=pv->idx;idxn=pvnext->idx;idxp=pvprev->idx;

		dx=ptptdistx(y[2*idxn],y[2*idxp]);
		dy=ptptdisty(y[2*idxn+1],y[2*idxp+1]);
		//dx=distx(pvnext,pvprev);
		//dy=disty(pvnext,pvprev);

		////// dH_dA * d2A_dr2 terms //////
		m->data[2*idx*m->tda+(2*idxn+1)]+=rho*0.5;
		m->data[2*idx*m->tda+(2*idxp+1)]-=rho*0.5;
		m->data[(2*idx+1)*m->tda+2*idxn]-=rho*0.5;
		m->data[(2*idx+1)*m->tda+2*idxp]+=rho*0.5;


		///// d2H_dA2 (dA_dr)^2 terms /////

		//I now think this is extra... j=i is taken care of in the loop
		//m->data[2*idx*m->tda+2*idx]+=0.5*d2E*dy*dy;
		//m->data[2*idx*m->tda+2*idx+1]-=0.5*d2E*dy*dx;
		//m->data[(2*idx+1)*m->tda+2*idx]-=0.5*d2E*dx*dy;
		//m->data[(2*idx+1)*m->tda+2*idx+1]+=0.5*d2E*dx*dx;

		for(j=0;j<nbneighb;j++){
			pvo=pcell->vertexlist[j];
			pvoprev=pcell->vertexlist[(j+nbneighb-1)%nbneighb];
			pvonext=pcell->vertexlist[(j+1)%nbneighb];
			idxo=pvo->idx;idxon=pvonext->idx;idxop=pvoprev->idx;

			dxo=ptptdistx(y[2*idxon],y[2*idxop]);
			dyo=ptptdisty(y[2*idxon+1],y[2*idxop+1]);
			//dxo=distx(pvonext,pvoprev);
			//dyo=disty(pvonext,pvoprev);

			m->data[2*idx*m->tda+2*idxo]+=0.5*d2E*dy*dyo;
			m->data[2*idx*m->tda+(2*idxo+1)]-=0.5*d2E*dy*dxo;
			m->data[(2*idx+1)*m->tda+2*idxo]-=0.5*d2E*dx*dyo;
			m->data[(2*idx+1)*m->tda+(2*idxo+1)]+=0.5*d2E*dx*dxo;
		}
	}
}

int hessian_area3d(gsl_matrix *m,int cellindex,const double y[]){

	int a,nbneighb,i,j;
	int idx,idxn,idxp,idxo,idxop,idxon;
	double dx,dy,dz,dxo,dyo,dzo;
	double rho,d2E;
	Dual *pcell;
	Node *pv,*pvprev,*pvnext;
	Node *pvo,*pvoprev,*pvonext;

	a = cellindex;pcell=web_dual[a];
	rho = pressure(pcell);
	d2E = -1.0*d2E_dA2(pcell);
	nbneighb = pcell->nb_vertices;

	for(i=0;i<nbneighb;i++){
		pv=pcell->vertexlist[i];
		pvprev=pcell->vertexlist[(i+1)%nbneighb];
		pvnext=pcell->vertexlist[(i+nbneighb-1)%nbneighb];
		idx=pv->idx;idxn=pvnext->idx;idxp=pvprev->idx;

		dx=ptptdistx(y[3*idxn],y[3*idxp]);
		dy=ptptdisty(y[3*idxn+1],y[3*idxp+1]);
		dz=ptptdistz(y[3*idxn+2],y[3*idxp+2]);
		//dx=distx(pvnext,pvprev);
		//dy=disty(pvnext,pvprev);
		//dz=distz(pvnext,pvprev);

		////// dH_dA * d2A_dr2 terms //////
		m->data[3*idx*m->tda+(3*idxn+1)]+=rho*0.5;
		m->data[3*idx*m->tda+(3*idxp+1)]-=rho*0.5;
		m->data[3*idx*m->tda+(3*idxn+2)]-=rho*0.5;
		m->data[3*idx*m->tda+(3*idxp+2)]+=rho*0.5;
		m->data[(3*idx+1)*m->tda+3*idxn]-=rho*0.5;
		m->data[(3*idx+1)*m->tda+3*idxp]+=rho*0.5;
		m->data[(3*idx+1)*m->tda+(3*idxn+2)]+=rho*0.5;
		m->data[(3*idx+1)*m->tda+(3*idxp+2)]-=rho*0.5;
		m->data[(3*idx+2)*m->tda+3*idxn]+=rho*0.5;
		m->data[(3*idx+2)*m->tda+3*idxp]-=rho*0.5;
		m->data[(3*idx+2)*m->tda+(3*idxn+1)]-=rho*0.5;
		m->data[(3*idx+2)*m->tda+(3*idxp+1)]+=rho*0.5;


		///// d2H_dA2 (dA_dr)^2 terms /////
		for(j=0;j<nbneighb;j++){
			pvo=pcell->vertexlist[j];
			pvoprev=pcell->vertexlist[(j+nbneighb-1)%nbneighb];
			pvonext=pcell->vertexlist[(j+1)%nbneighb];
			idxo=pvo->idx;idxop=pvoprev->idx;idxon=pvonext->idx;

			dxo=ptptdistx(y[3*idxon],y[3*idxop]);
			dyo=ptptdisty(y[3*idxon+1],y[3*idxop+1]);
			dzo=ptptdistz(y[3*idxon+2],y[3*idxop+2]);
			//dxo=distx(pvonext,pvoprev);
			//dyo=disty(pvonext,pvoprev);
			//dzo=distz(pvonext,pvoprev);

			m->data[3*idx*m->tda+3*idxo]+=0.5*d2E*dy*dyo;
			m->data[3*idx*m->tda+(3*idxo+1)]-=0.5*d2E*dy*dxo;
			m->data[(3*idx+1)*m->tda+3*idxo]-=0.5*d2E*dx*dyo;
			m->data[(3*idx+1)*m->tda+(3*idxo+1)]+=0.5*d2E*dx*dxo;
		}
	}
}

int test_symmetric(gsl_matrix *m){

	double diff;
	int i,j;

	if(m->size1!=m->size2){printf("Matrix is not square!\n");return 0;}
	for(i=0;i<m->size1;i++){
		for(j=0;j<m->size2;j++){
			diff=gsl_matrix_get(m,i,j)-gsl_matrix_get(m,j,i);
			if(abs(diff)>0.0000001){
				printf("Matrix (i,j) entry is not symmetric with diff=%lf",diff);fflush(stdout);
			}
		}
	}
	return 1;
}

int test_eig(gsl_matrix *M,gsl_vector *eval,gsl_matrix *evec){

	if(M->size1 != M->size2){locerror("test_eig()","matrix is not square!");}
	if(eval->size!=M->size1){locerror("test_eig()","output vector not correct size!");}
	if(evec->size1!=M->size1 || evec->size2!=M->size2){locerror("test_eig()","output and input matrix are differently sized!");}
	
	gsl_matrix *m = gsl_matrix_alloc(M->size1,M->size2);
	gsl_matrix_memcpy(m,M);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(m->size1);

	gsl_eigen_symmv(m,eval,evec,w);

	gsl_eigen_symmv_free(w);gsl_matrix_free(m);

	gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);
}

void evol(const double delta_t,const char *method){
    int i,count,filecount;
    FILE *ofp;
    char filename[50];
	double *y;
	int numdim=NUMDIM;
	double p[1]={(double) MUP};

	if(strcmp(method,"all contract")==0){}
	else if(strcmp(method,"mechanical contract")==0){}
	else if(strcmp(method,"")==0){}
	else{}

	y = (double *)malloc(numdim*nb_vertex_tot*sizeof(double));
	for(i=0;i<nb_vertex_tot;i++){
		y[numdim*i]=web[i]->x;
		y[numdim*i+1]=web[i]->y;
		if(numdim==3){y[numdim*i+2]=web[i]->z;}
	}

    
    //p=allocevolparam();
    const gsl_odeiv2_step_type *type_bsimp=gsl_odeiv2_step_bsimp; //Implicit Burlish-Stoer method
	const gsl_odeiv2_step_type *type_rkf45=gsl_odeiv2_step_rkf45; //Explicit embedded RK-Fehlberg (4,5) method
    gsl_odeiv2_step *s=gsl_odeiv2_step_alloc(type_rkf45,numdim*nb_vertex_tot);
    gsl_odeiv2_control *c=gsl_odeiv2_control_y_new(1e-5,0.0);
    gsl_odeiv2_evolve *e=gsl_odeiv2_evolve_alloc(numdim*nb_vertex_tot);

    gsl_odeiv2_system sys = {evolfunc,evoljac,numdim*nb_vertex_tot,p};

    double t=0.0,t1=delta_t,t0=0.0;
    double h=1e-3;

	//printf("starting iterations....\n");fflush(stdout);
    count=0;filecount=0;
    while(t<t1){
		t0=t;
        int status=gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,t1,&h,y);
		check_cell_state_stress_activation(25.0,2.0,t-t0); //toroidal
		//check_cell_state_stress_activation(4.0,1.0,t-t0); //spherical
        if(status!=GSL_SUCCESS){break;printf("failure :(!\n");fflush(stdout);}
        count++;
		if(count%1000==0){
			printf("time = %lf. ",t);fflush(stdout);
            //for(i=0;i<nb_cell_tot;i++){printf("%.4e %.4e %.4e\n",y[3*i],y[3*i+1],y[3*i+2]);fflush(stdout);}
		}

    }
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

	for(i=0;i<nb_vertex_tot;i++){
		web[i]->x=y[numdim*i];
		web[i]->y=y[numdim*i+1];
		if(numdim==3){web[i]->z=y[numdim*i+2];}
		else{web[i]->z=0.0;}
	}
	free(y);
}

void evol_stress_activation(const double delta_t){
    int i,count,filecount;
    FILE *ofp;
    char filename[50];
	double *y;
	int numdim=NUMDIM;
	double p[1]={(double) MUP};

	y = (double *)malloc(numdim*nb_vertex_tot*sizeof(double));
	for(i=0;i<nb_vertex_tot;i++){
		y[numdim*i]=web[i]->x;
		y[numdim*i+1]=web[i]->y;
		if(numdim==3){y[numdim*i+2]=web[i]->z;}
	}

    
    //p=allocevolparam();
    const gsl_odeiv2_step_type *type_bsimp=gsl_odeiv2_step_bsimp; //Implicit Burlish-Stoer method
	const gsl_odeiv2_step_type *type_rkf45=gsl_odeiv2_step_rkf45; //Explicit embedded RK-Fehlberg (4,5) method
    gsl_odeiv2_step *s=gsl_odeiv2_step_alloc(type_rkf45,numdim*nb_vertex_tot);
    gsl_odeiv2_control *c=gsl_odeiv2_control_y_new(1e-5,0.0);
    gsl_odeiv2_evolve *e=gsl_odeiv2_evolve_alloc(numdim*nb_vertex_tot);

    gsl_odeiv2_system sys = {evolfunc,evoljac,numdim*nb_vertex_tot,p};

    double t=0.0,t1=delta_t,t0=0.0;
    double h=1e-3;

	//printf("starting iterations....\n");fflush(stdout);
    count=0;filecount=0;
    while(t<t1){
		t0=t;
        int status=gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,t1,&h,y);
		check_cell_state_stress_activation(25.0,2.0,t-t0); //toroidal
		//check_cell_state_stress_activation(4.0,1.0,t-t0); //spherical
        if(status!=GSL_SUCCESS){break;printf("failure :(!\n");fflush(stdout);}
        count++;
		if(count%1000==0){
			printf("time = %lf. ",t);fflush(stdout);
            //for(i=0;i<nb_cell_tot;i++){printf("%.4e %.4e %.4e\n",y[3*i],y[3*i+1],y[3*i+2]);fflush(stdout);}
		}

    }
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

	for(i=0;i<nb_vertex_tot;i++){
		web[i]->x=y[numdim*i];
		web[i]->y=y[numdim*i+1];
		if(numdim==3){web[i]->z=y[numdim*i+2];}
		else{web[i]->z=0.0;}
	}
	free(y);
}

void lead_contract(int *list,int length){

	int i,start=0;


	for(i=0;i<length;i++){
		if(TOROIDAL==0 && list[i]==0){printf("oops, should not modify cell0 for non-toroidal!\n");fflush(stdout);}
		printf("modifying cell %d!\n",list[i]);fflush(stdout);
		web_dual[list[i]]->time=0.2; //time in contracting state
		web_dual[list[i]]->marker.clone_index=1; //set contracting state
		web_dual[list[i]]->area_soll=0.4*web_dual[list[i]]->area_soll;
	}
}

void init_evol_stress_activation(Dual *pcell){

	int i;

	for(i=0;i<nb_bond_tot;i++){
		web_bond[i]->kappa=0.;
		web_bond[i]->L0=1.0;
		web_bond[i]->lambda=0.0;
	}

	int list[3]={1,2,3};
	lead_contract(list,3);

}

void check_cell_state_stress_activation(double base,double thresh,double dt){

	int i,start=0;
	double str_det,str_tr;
	matrix22 str;
	Dual *pcell;

	if(TOROIDAL==0){start=1;}

	for(i=start;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		str=stress_cells_2d(&i,1);
		str_tr=str.c[0][0]+str.c[1][1];
		str_det=str.c[0][0]*str.c[1][1]-str.c[0][1]*str.c[1][0];

		pcell->time-=dt; //update cell clock

		if(abs(str_det-base)>thresh && pcell->marker.clone_index==0){
			pcell->time=0.2; //time in contracting state
			pcell->marker.clone_index=1; //set contracting state
			pcell->area_soll = 0.4*pcell->area_soll; //modify energetics
		} //stress condition

		if(pcell->marker.clone_index==1 && pcell->time<0.0){
			pcell->time = 1.0; //time in refractory state
			pcell->marker.clone_index=2; //set refractory state
			pcell->area_soll = pcell->area_soll0; //modify energetics
		}

		if(pcell->marker.clone_index==2 && pcell->time<0.0){
			pcell->time = 0.; //set t=0, shouldn't matter
			pcell->marker.clone_index=0; //set primed state
		}

	}

}

void init_evol_total_activation(Dual *pcell){

	int i;
	Bond *pbond;

	//do nothing with bond energies yet
	for(i=0;i<nb_bond_tot;i++){
		web_bond[i]->kappa=0.;
		web_bond[i]->L0=0.;
		web_bond[i]->lambda=0.;
	}
	
	pcell->marker.border=1;

	for(i=0;i<nb_cell_tot;i++){
		web_dual[i]->area_soll=0.5;
		web_dual[i]->kA=6.0;
		web_dual[i]->sigma=20.;
	}

}


void evol_total_activation(const double delta_t,const double E0){
    int i,count,filecount;
    FILE *ofp;
    char filename[50];
	double *y;
	int numdim=NUMDIM;
	double p[1]={(double) E0};

	y = (double *)malloc(numdim*nb_vertex_tot*sizeof(double));
	for(i=0;i<nb_vertex_tot;i++){
		y[numdim*i]=web[i]->x;
		y[numdim*i+1]=web[i]->y;
		if(numdim==3){y[numdim*i+2]=web[i]->z;}
	}

    //p=allocevolparam();
    const gsl_odeiv2_step_type *type_bsimp=gsl_odeiv2_step_bsimp; //Implicit Burlish-Stoer method
	const gsl_odeiv2_step_type *type_rkf45=gsl_odeiv2_step_rkf45; //Explicit embedded RK-Fehlberg (4,5) method
    gsl_odeiv2_step *s=gsl_odeiv2_step_alloc(type_rkf45,numdim*nb_vertex_tot);
    gsl_odeiv2_control *c=gsl_odeiv2_control_y_new(1e-5,0.0);
    gsl_odeiv2_evolve *e=gsl_odeiv2_evolve_alloc(numdim*nb_vertex_tot);

    gsl_odeiv2_system sys = {evolfunc,evoljac,numdim*nb_vertex_tot,p};

    double t=0.0,t1=delta_t,t0=0.0;
    double h=1e-3;

	//printf("starting iterations....\n");fflush(stdout);
    count=0;filecount=0;
    while(t<t1){
		t0=t;
        int status=gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,t1,&h,y);
		//no need to change cell states
		//check_cell_state(25.0,2.0,t-t0); //toroidal
		//check_cell_state(4.0,1.0,t-t0); //spherical
        if(status!=GSL_SUCCESS){break;printf("failure :(!\n");fflush(stdout);}
        count++;
		if(count%1000==0){
			printf("time = %lf. ",t);fflush(stdout);
            //for(i=0;i<nb_cell_tot;i++){printf("%.4e %.4e %.4e\n",y[3*i],y[3*i+1],y[3*i+2]);fflush(stdout);}
		}

    }
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

	for(i=0;i<nb_vertex_tot;i++){
		web[i]->x=y[numdim*i];
		web[i]->y=y[numdim*i+1];
		if(numdim==3){web[i]->z=y[numdim*i+2];}
		else{web[i]->z=0.0;}
	}
	free(y);
}

void posfluc(double sigma){
	int i,j;
	double *temp;
	double dr;
	Node *pv;

	if(sigma>1e-6){
		temp=(double*)calloc(nb_vertex_tot,sizeof(double));
		for(i=0;i<nb_vertex_tot;i++){
			pv=web[i];
			for(j=0;j<3;j++){temp[i]+=dist(pv->pneighb[j],web[i])/3.;}
			temp[i]=0.8; //hard coded in for now, should change
		}
		for(i=0;i<nb_vertex_tot;i++){
			dr=gsl_ran_gaussian(rng,sigma);
			web[i]->x+=dr*temp[i];
			dr=gsl_ran_gaussian(rng,sigma);
			web[i]->y+=dr*temp[i];
			if(NUMDIM==3){dr=gsl_ran_gaussian(rng,sigma);web[i]->z+=dr*temp[i];}
		}
		free(temp);
	}
}
