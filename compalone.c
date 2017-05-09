/* Standalone program (does not need opengl) */
// export CPPFLAGS="-I/sw/include/"
// export LDFLAGS="-L/sw/lib"

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>

#include "const.h"
#include "lattice.h"
#include "save.h"
#include "read.h"
#include "pngwrite.h"
#include "locerror.h"
#include "measurements.h"
#include "evolution.h"
#include "compalone.h"
#include "opti.h"
#include "in_development.h"
#include "dynamics.h"


char outfiledir[200];

//prototype the readcommandline() function, since it's quite long
void readcommandline(int argc,char **argv);
//prototype the checkOutputDir() function
void checkOutputDir();

int divc_selective(){
	int number=0,number2=0,nbvert=0;
	int i,j,k,retval=0,rand2;
	double rand,bound,mult_factor=0.80,tot;//variables for clone cells on the boundary
	double bound2,factor2=1.2;//variables for clone cells inside boundary
	Dual *pcell;

	for(i=1;i<nb_cell_tot;i++){
		pcell = web_dual[i];
		if(pcell->marker.clone_index==1){number++;}
	}
	rand = (gsl_rng_uniform(rng));
	//if(number==0){mult_factor=0.0;/*division weight for boundary clones*/}
	//if(number2==0){factor2=0.0;/*division weight for non-boundary clones*/}
	tot = (mult_factor-1.0)*number+(factor2-1.0)*number2+nb_cell_tot; //total weight
	bound = mult_factor*number/tot;
	bound2 = factor2*number2/tot + bound;
	//following if statements test for which type of cell is dividing and then picks
	//a random one of those cells at a random vertex and calls division_cell_slow()
	if(rand<bound){
		while(nbvert<=4){
			rand2 = ceil((gsl_rng_uniform(rng)*(double)number));
			j=0;
			for(i=0;i<nb_cell_tot;i++){
				pcell = web_dual[i];
				//if(pcell->marker.border==1){j++;}
				if(pcell->marker.clone_index==1){j++;}
				if(j==rand2){
					nbvert=pcell->nb_vertices;
					if(nbvert>4){division_cell_slow(pcell,floor(gsl_rng_uniform(rng)*nbvert));retval=1;break;}
				}
			}
		}
	}
	if(rand>bound && rand<bound2){
		while(nbvert<=4){
			rand2 = ceil((gsl_rng_uniform(rng)*(double)(number2)));
			j=0;
			for(i=0;i<nb_cell_tot;i++){
				pcell = web_dual[i];
				if(pcell->marker.clone_index==0){j++;}
				if(j==rand2){
					nbvert=pcell->nb_vertices;
					if(nbvert>4){division_cell_slow(pcell,floor(gsl_rng_uniform(rng)*nbvert));retval=1;break;}
				}
			}
		}
	}
	if(rand>bound2){
		while(nbvert<=4){
			rand2 = ceil((gsl_rng_uniform(rng)*(double)(nb_cell_tot-number-number2-1)));
			j=0;
			for(i=1;i<nb_cell_tot;i++){
				pcell = web_dual[i];
				if(pcell->marker.border==0 && pcell->marker.tension_index==1.0){j++;}
				if(j==rand2){
					nbvert=pcell->nb_vertices;
					if(nbvert>4){division_cell_slow(pcell,floor(gsl_rng_uniform(rng)*nbvert));retval=1;break;}
				}
			}
		}
	}
	if(retval!=1){printf("number = %d, number2 = %d, rand2 = %d, rand = %lf, bound = %lf, bound2 = %lf, nb_cell_tot = %d\n",number,number2,rand2,rand,bound,bound2,nb_cell_tot);}//in case of failure for diagnostics
	return retval;
}

int divc_stocha()
{
	int i;
	int retval=0;
	Dual *pcell;


	i = 1 + (int) (gsl_rng_uniform(rng)*(nb_cell_tot-1));
	pcell=web_dual[i];
	division_cell_slow(pcell,(int)(gsl_rng_uniform(rng)*pcell->nb_vertices));
	retval=1;
	return retval;
}

unsigned int nb_killed=0;

void heartevol()
{
	int i, j, k;
	int nb_lat=1;
	char filename[256],svgfname[256],infofname[256];
	char svgdir[30],infodir[30];
	FILE *pfile1, *pfile2;
	int numtimesteps=1000;
	Dual *pcell_mid;

	for (i=0;i<nb_lat;i++){
		checkOutputDir(); //make sure our output directory is sensible
		sprintf(svgdir,"%s/svgs",outfiledir);sprintf(infodir,"%s/info",outfiledir);
        //init_lattice();
		init_lattice_torus();
		updatecellbondneighbs();
        printf("we have %d cells to play with!\n",nb_cell_tot);

		sprintf(filename,"initial");
		createFilenames(256,svgfname,infofname,svgdir,infodir,filename);
		svglattice(svgfname,2);infoprint(infofname,2);

		optimize(); //relax 
		//set parameters of this center cell to have no energetic contribution
		pcell_mid=findCenterCell();
		pcell_mid->marker.border=1;//this is a "border cell"
		//to be excluded from energetic considerations
		optimize(); //put in relaxed state after changing middle cell

		//let's check out the relaxed state
		sprintf(filename,"initial_relaxed");
		createFilenames(256,svgfname,infofname,svgdir,infodir,filename);
		svglattice(svgfname,2);infoprint(infofname,2);

		//now change energetics so everything is more tensile
		init_evol_total_activation(pcell_mid); 
		//initialize thermal fluction parameters
		double E0=0.;
		//and then evolve and output results
		for(j=0;j<numtimesteps;j++){
			evol_total_activation(0.001,E0);
			printf("%d%%\r",(int)(100*j/numtimesteps));fflush(stdout);
			sprintf(filename,"%.4i",j);
			createFilenames(256,svgfname,infofname,svgdir,infodir,filename);
			svglattice(svgfname,2);infoprint(infofname,2);
			//if(j%300==0){lead_contract(list,3);}
		}
		printf("\n");

		for(k=0;k<3;k++){free(forcevector[k]);}
		free(forcevector);
		forcevector=(double**)malloc(3*sizeof(double*));
		for(k=0;k<3;k++){
			forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
		}

	}

	kill_lattice();
	for(k=0;k<3;k++){free(forcevector[k]);}
	free(forcevector);
}

GlobalParams *setGlobalModelParameters(){
	int i;
	GlobalParams *p;

	//Rest length parameters
	p->A0=LX*LY/((double) nb_cell_tot);
	p->P0=3.73*p->A0;
	p->L0=p->P0/6.;
	//Prefactors for quadratic models
	p->kA=1.;
	p->kP=1.;
	p->kL=0.;
	//Prefactors for linear models
	p->sigma=0.;
	p->lambda=0.;

	for(i=0;i<nb_cell_tot;i++){
		web_dual[i]->area_soll=p->A0;
		web_dual[i]->area_soll0=p->A0;
		web_dual[i]->perim_soll=p->P0;
		web_dual[i]->sigma = p->sigma;
		web_dual[i]->kP = p->kP;
		web_dual[i]->kA = p->kA;
	}
	for(i=0;i<nb_bond_tot;i++){
		web_bond[i]->kappa=p->kL;
		web_bond[i]->lambda=p->lambda;
		web_bond[i]->L0=p->L0;
	}

}

void rigidityevol(){
	GlobalParams *p;
	char filename[256],svgfname[256],infofname[256];
	char svgdir[50],infodir[50];

	checkOutputDir(); //make sure our output directory is sensible
	sprintf(svgdir,"%s/svgs",outfiledir);sprintf(infodir,"%s/info",outfiledir);
	createFilenames(256,svgfname,infofname,svgdir,infodir,filename);

	init_lattice_torus(); //creates lattice and basic data structures
	updatecellbondneighbs(); //creates and organizes bond structure
	p=setGlobalModelParameters();
	optimize(); //updates geometry information and finds energetic minimum


	kill_lattice();
	for(int k=0;k<3;k++){free(forcevector[k]);}
	free(forcevector);
}

gsl_rng *rng;

// Main function, primarily to call content functions and set certain global variables
int main(int argc, char *argv[]){

	//initialize model type here, potentially redefined by input
	strcpy(AREAMODELTYPE,"quadratic"); //options: "none", "quadratic", and "hyperelastic" for now
	strcpy(PERIMMODELTYPE,"quadratic"); //options: "none", "linear", and "quadratic" for now
	//initialize output directory
	strcpy(outfiledir,"./areamodel"); 

	readcommandline(argc,argv);
	checkModelTypes();

	rng=gsl_rng_alloc(gsl_rng_default);
	/* The following functions can be activated in compalone.h */
	rigidityevol();
	return 0;
}

void readcommandline(int argc,char *argv[]){

	int i;
	char temp[200];

	if(argc>1){
		for(i=1;i<argc;i++){
			if(!strcmp(argv[i],"-f")){
				if(i+1==argc){locerror("readcommandline()","No filename!");}
				i++;snprintf(temp,sizeof(temp),"./%s",argv[i]);
				strncpy(outfiledir,temp,sizeof(outfiledir));
			}
			else if(!strcmp(argv[i],"-areamodeltype")){
				if(i+1==argc){locerror("readcommandline()","No area model type!");}
				i++;snprintf(temp,sizeof(temp),"./%s",argv[i]);
				strncpy(AREAMODELTYPE,temp,sizeof(AREAMODELTYPE));
			}
			else if(!strcmp(argv[i],"-perimmodeltype")){
				if(i+1==argc){locerror("readcommandline()","No perimeter model type!");}
				i++;strncpy(PERIMMODELTYPE,argv[i],sizeof(PERIMMODELTYPE));
			}


		}
	}
	else{
		//could do initializations here as well
		//didn't because of code organization I guess

	}
}

void checkOutputDir(){
	int retval;
	char infodir[205],svgdir[205],jpgdir[205];
	struct stat st={0};
	//utilize system calls to do directory stuff, 
	snprintf(infodir,sizeof(infodir),"%s/info",outfiledir);
	snprintf(svgdir,sizeof(svgdir),"%s/info",outfiledir);
	snprintf(jpgdir,sizeof(jpgdir),"%s/info",outfiledir);
	
	//make top level data directory if it doesn't exist
	if(stat(outfiledir,&st)==-1) retval=mkdir(outfiledir,0755);
	//make other data directories if they don't exist
	if(stat(infodir, &st)==-1) retval=mkdir(infodir,0755);
	if(stat(jpgdir,&st)==-1) retval=mkdir(jpgdir,0755);
	if(stat(svgdir,&st)==-1) retval=mkdir(svgdir,0755);
}
