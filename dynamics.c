#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "definelat.h"
#include "locerror.h"
#include "lattice.h"
#include "measurements.h"
#include "dynamics.h"
#include "compalone.h"

Lat **web_lat;

Parameter *allocparam(){
    Parameter *a;
    a=(Parameter *)malloc(sizeof(Parameter));
    return a;
}
Lat **createlatvector(int dim){
	Lat **a;
	a=(Lat**)malloc(dim*sizeof(Lat*));
	if(!a){locerror("createlatvector","out of memory!");}
	return a;
}

Lat **realloclatvector(Lat **latvec,int dim){
	Lat **a;
	a = (Lat**)realloc(latvec,dim*sizeof(Lat*));
	if(!a) locerror("realloclatvector","out of memory!");
	return a;
}

Lat *alloclat(){
	Lat *a;
	a=(Lat*)malloc(sizeof(Lat));
	if(!a){locerror("alloclat","out of memory!");}
	return a;
}
int averageND(const Dual *pcell,const double y[],double *avgN, double *avgD){
    int i,j,nbneighb;

    *avgN=0.0;*avgD=0.0;
    nbneighb=pcell->nb_vertices;
    for(i=0;i<nbneighb;i++){
        j=indexcellinwebdual(pcell->celllist[i]);
        *avgN+=y[3*j]/nbneighb;
        *avgD+=y[3*j+1]/nbneighb;
    }
    return 0;
}

int func(double t, const double y[], double f[], void *par){
    int i,j,k,nbneighb;
    Dual *pcell;
    double avgN,avgD;
    double N,D,R;
    double BN,BD,T,KC,KR,BR;
    double *result;
    Parameter *params = (Parameter *)par;

    BN=params->bn;
    BD=params->bd;
    T=params->t;
    KC=params->kc;
    KR=params->kr;
    BR=params->br;
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        nbneighb=pcell->nb_vertices;
        N=y[3*i];D=y[3*i+1];R=y[3*i+2];
        //find the terms that depend on neighbors
        averageND(pcell,y,&avgN,&avgD);
        f[3*i]=1/T*(BN - N*(1+avgD+D/KC));
        f[3*i+1]=1/T*(BD/(1+R) - D*(1+avgN +N/KC));
        f[3*i+2]= BR*(gsl_pow_3(N*avgD)/(gsl_pow_3(KR) + gsl_pow_3(N*avgD))) - R;
    }
    return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy, double dfdt[], void *par){
    int i,j,k,nbneighb;
    Dual *pcell;
    double avgN,avgD;
    double N,D,R;
    double BN,BD,T,KC,KR,BR;
    double denom,numer,sum;
    double *result;
    Parameter *params = (Parameter *)par;

    BN=params->bn;
    BD=params->bd;
    T=params->t;
    KC=params->kc;
    KR=params->kr;
    BR=params->br;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,3*nb_cell_tot,3*nb_cell_tot);
    gsl_matrix *m = &dfdy_mat.matrix;

    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        nbneighb=pcell->nb_vertices;

        N=y[3*i];D=y[3*i+1];R=y[3*i+2];
        //find the terms that depend on neighbors
        averageND(pcell,y,&avgN,&avgD);
        denom=gsl_pow_3(KR)+gsl_pow_3(N*avgD);
        numer=gsl_pow_2(N*avgD);
        sum=(numer*N/denom-gsl_pow_3(N*avgD)*numer*N/gsl_pow_2(denom))/nbneighb;
        //trans terms in jacobian
        for(k=0;k<nbneighb;k++){
            j=indexcellinwebdual(pcell->celllist[k]);
            gsl_matrix_set(m,3*i   ,3*j+1 ,-N/nbneighb/T);
            gsl_matrix_set(m,3*i+1 ,3*j   ,-D/nbneighb/T);
            gsl_matrix_set(m,3*i+2 ,3*j+1 ,BR*sum);
        }
        //cis terms in jacobian
        gsl_matrix_set(m,3*i   ,3*i   ,-(1+avgD+D/KC)/T);
        gsl_matrix_set(m,3*i   ,3*i+1 ,-N/KC/T);
        gsl_matrix_set(m,3*i+1 ,3*i   ,-D/KC/T);
        gsl_matrix_set(m,3*i+1 ,3*i+1 ,-(1+avgN+N/KC)/T);
        gsl_matrix_set(m,3*i+1 ,3*i+2 ,-BD/(gsl_pow_2(1+R)*T));
        sum=numer*avgD/denom-gsl_pow_3(N*avgD)*numer*avgD/gsl_pow_2(denom);
        gsl_matrix_set(m,3*i+2 ,3*i   ,BR*sum);
        gsl_matrix_set(m,3*i+2 ,3*i+2 ,-1);
        dfdt[i]=0.0;
    }
    return GSL_SUCCESS;
}

void dynamics(){
    int i,count,filecount;
    Parameter *p;
    FILE *ofp;
    char filename[50];
	double *y;
	Lat *plat;

	y = (double *)malloc(3*nb_cell_tot*sizeof(double));
	for(i=0;i<nb_cell_tot;i++){
		plat=web_lat[i];
		y[3*i]=web_lat[i]->notch;
		y[3*i+1]=web_lat[i]->delta;
		y[3*i+2]=web_lat[i]->reg;
	}
    
    p=allocparam();
    const gsl_odeiv_step_type *type=gsl_odeiv_step_bsimp;
    gsl_odeiv_step *s=gsl_odeiv_step_alloc(type,3*nb_cell_tot);
    gsl_odeiv_control *c=gsl_odeiv_control_y_new(1e-6,0.0);
    gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc(3*nb_cell_tot);

    p->bn=102.0;p->bd=98.0;p->t=1.0;p->kc=0.5;p->kr=10.0;p->br=100.0;
    gsl_odeiv_system sys = {func,jac,3*nb_cell_tot,p};

    double t=0.0,t1=4.0;
    double h=1e-6;

	//printf("starting iterations....\n");fflush(stdout);
    count=0;filecount=0;
    while(t<t1){
        int status=gsl_odeiv_evolve_apply(e,c,s,&sys,&t,t1,&h,y);
        if(status!=GSL_SUCCESS){break;printf("failure :(!\n");fflush(stdout);}
/*
        if((count%4)==0){ //prints diagnostic files
            sprintf(filename,"txt/dynamics_%.4d.txt",filecount);
            ofp=fopen(filename,"w");
            fprintf(ofp,"%lf\n",t);
            for(i=0;i<nb_cell_tot;i++){fprintf(ofp,"%.4e %.4e %.4e\n",y[3*i],y[3*i+1],y[3*i+2]);}
            fclose(ofp);
            printf("printed file num %d!\n",filecount);fflush(stdout);
            filecount++;
        }
*/

        count++;
		if(count%1000==0){
			printf("time = %lf. ",t);fflush(stdout);
            for(i=0;i<nb_cell_tot;i++){printf("%.4e %.4e %.4e\n",y[3*i],y[3*i+1],y[3*i+2]);fflush(stdout);}
		}

    }
	printf("\n");fflush(stdout);
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);

	for(i=0;i<nb_cell_tot;i++){
		web_lat[i]->notch=y[3*i];
		web_lat[i]->delta=y[3*i+1];
		web_lat[i]->reg=y[3*i+2];
	}
	free(y);
}

void init_dynamics(){
	Lat *plat;
	int i;

	web_lat = createlatvector(nb_cell_tot);
	for(i=0;i<nb_cell_tot;i++){web_lat[i] = alloclat();}
	
	for(i=0;i<nb_cell_tot;i++){
		plat=web_lat[i];

		plat->notch=0.5;
		plat->delta=1.0;
		plat->reg=0.3;
	}
	web_lat[2]->notch=0.1;
	web_lat[2]->delta=10.5;
	web_lat[2]->reg=1.0;

}

void init_dynamics_rnd(int init){
	Lat *plat;
	int i;

	if(init==1){
		web_lat = createlatvector(nb_cell_tot);
		for(i=0;i<nb_cell_tot;i++){web_lat[i] = alloclat();}
	}

	for(i=0;i<nb_cell_tot;i++){
		plat=web_lat[i];

		plat->notch=0.5;
		plat->delta=1.0;
		plat->reg=0.3;
	}
	i = floor(nb_cell_tot*gsl_rng_uniform(rng));
	web_lat[i]->notch=0.1;
	web_lat[i]->delta=10.5;
	web_lat[i]->reg=1.0;
}

void setclone_dynamics(){

	int i;
	double avgN=0.,avgD=0.;
	Dual *pcell;

	for(i=0;i<nb_cell_tot;i++){
		avgN+=(web_lat[i])->notch/nb_cell_tot;
		avgD+=(web_lat[i])->delta/nb_cell_tot;
	}

	//printf("avgN = %lf, avgD = %lf\n",avgN,avgD);
	for(i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		if((web_lat[i])->notch > avgN && (web_lat[i])->delta < avgD){pcell->marker.clone_index=0;}
		else if((web_lat[i])->notch < avgN && (web_lat[i])->delta > avgD){pcell->marker.clone_index=1;}
	}
}

void division_dynamics(int origcell){

	int i;
	Lat **web_lat_temp;

	web_lat_temp = createlatvector(nb_cell_tot-1);
	for(i=0;i<nb_cell_tot-1;i++){web_lat_temp[i] = alloclat();}
	for(i=0;i<nb_cell_tot-1;i++){
		web_lat_temp[i]->notch = web_lat[i]->notch;
		web_lat_temp[i]->delta = web_lat[i]->delta;
		web_lat_temp[i]->reg = web_lat[i]->reg;
	}

	web_lat = realloclatvector(web_lat,nb_cell_tot);
	web_lat[nb_cell_tot-1]=alloclat();

	for(i=0;i<nb_cell_tot-1;i++){
		web_lat[i]->notch = web_lat_temp[i]->notch;
		web_lat[i]->delta = web_lat_temp[i]->delta;
		web_lat[i]->reg = web_lat_temp[i]->reg;
	}

	web_lat[nb_cell_tot-1]->notch = (web_lat_temp[origcell]->notch)*1.1;
	web_lat[nb_cell_tot-1]->delta = (web_lat_temp[origcell]->delta)*0.9;
	web_lat[nb_cell_tot-1]->reg = (web_lat_temp[origcell]->reg)*0.9;

	for(i=0;i<nb_cell_tot-1;i++){free(web_lat_temp[i]);}
	free(web_lat_temp);
}

int hstate_check(){

	int i,j,start,status;
	Dual *pcell,*pncell;

	if(TOROIDAL>0){start=0;}
	else{start=1;}

	for(i=start;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		if(pcell->marker.clone_index==1){
			for(j=0;j<pcell->nb_vertices;j++){
				pncell=pcell->celllist[j];
				if(pncell->marker.clone_index==1){pcell->marker.clone_index=0;}
			}
		}
		else if(pcell->marker.clone_index==0){
			status=0;
			for(j=0;j<pcell->nb_vertices;j++){
				pncell=pcell->celllist[j];
				if(pncell->marker.clone_index==1){status+=1;}
			}
			if(status==0){pcell->marker.clone_index=1;}
		}
		else{printf("We have clone index = %d!\n",pcell->marker.clone_index);fflush(stdout);}
	}
}
