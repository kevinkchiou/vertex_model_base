#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multimin.h>

#include "const.h"
#include "definelat.h"
#include "lattice.h"
#include "measurements.h"
#include "read.h"
#include "save.h"
#include "compalone.h"
#include "locerror.h"
#include "in_development.h"
#include "opti.h"



//gives each vertex an offset
void coord_mod(){
    int i;
    Node *pvert;

    for(i=0;i<nb_vertex_tot;i++){
       pvert = web[i];
       pvert->z = 12.0;
    }
}

void border_check(){
   int i,j,k;
   Dual *pcell,*pcellnb;
   Node *pvert;

   for(i=0;i<nb_vertex_tot;i++){
      pvert = web[i];
      pvert->border = 0;
   }
   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      pcell->marker.border=0;
   }

   /*
   for(i=start;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      if((pcell->marker.tension_index)>1.001){
         pcell->area_soll=0.8;
         pcell->area_soll0=pcell->area_soll;
         for(j=0;j<(pcell->nb_vertices);j++){
            pvert=pcell->vertexlist[j];
            for(k=0;k<3;k++){
               pcellnb = pvert->pncell[k];
               if((pcellnb->marker.tension_index)<1.001){
                  pvert->border=1;
                  pcell->marker.border=1;
                  pcell->area_soll=0.8;
                  pcell->area_soll0=pcell->area_soll;
               }
            }
         }
      }
   }
   */
}

void division_cell_slow(Dual *pcell, int ivertex){
   int i,max=3,icell;//max is # of steps to get to max size
   double starting_area;
   Node *pvertex;

   pvertex=pcell->vertexlist[ivertex];
   icell=indexcell(pvertex,pcell);
   if(pcell->area_soll > pcell->area_soll0){starting_area = pcell->area_soll0;}
   else{starting_area = pcell->area_soll;}
   for(i=0;i<max;i++){
      pcell->area_soll=(1.0+((double) (i+1))/((double) max))*(starting_area);
      optimizet1();
   }
   division_eq(pvertex,icell);
}

void tension_modify(){
   Dual *pcell;
   int icell;

   printf("modifying tension...\n");
   add_factor = add_factor2;
   border_check();
   for(icell=0;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      if(pcell->marker.tension_index!=1.0){pcell->marker.tension_index+=(add_factor2-add_factor1);}
      border_check();
   }
   optimizet1();
}

void tension_modify1(){
	Dual *pcell;
	int icell,count=0,numdel;

	printf("modifying tension...\n");
	for(icell=0;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		if(pcell->marker.clone_index==1){pcell->marker.tension_index=8.0;pcell->area_soll=13.0;count++;}
		else{pcell->marker.tension_index=10.0;pcell->area_soll=9.0;}
		border_check();
	}
	printf("we have %d delta cells, and %d notch cells!\n",count,nb_cell_tot-count);
	optimize();
}

double tension_modify2(double eps){
   Dual *pcell;
   int icell;
   
   printf("modifying tension...\n");
   border_check();
   for(icell=0;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      pcell->marker.tension_index=1.0+2*eps*(gsl_rng_uniform(rng)-0.5);
   }
   optimizet1();
   return eps;
}

void tension_uniform(){
    int icell;
    Dual *pcell;

    printf("modifying tensions to attempt random lattice...");fflush(stdout);
    for(icell=0;icell<nb_cell_tot;icell++){
        pcell=web_dual[icell];
        pcell->area_soll=(LX*LY-1.0)/nb_cell_tot;
        pcell->marker.tension_index=pcell->area_soll;
    }
    optimizet1();
}

void tension_uniform_eraseclone(){
    int icell;
    Dual *pcell;

    printf("modifying tensions to attempt random lattice...");fflush(stdout);
    for(icell=0;icell<nb_cell_tot;icell++){
        pcell=web_dual[icell];
        pcell->marker.clone_index=0;
        pcell->area_soll=(LX*LY-1.0)/nb_cell_tot;
        pcell->marker.tension_index=pcell->area_soll;
    }
    optimizet1();
}

void kagome_pre_modify(){
	Dual *pcell,*pncell;
	int icell,jcell,deltaflag,nbneighb;

	printf("placing delta cells in cell with low numbers of neighbors...\n");
	for(icell=1;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		if(nbneighb<5){
			deltaflag=1;
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){deltaflag=0;}
			}
			if(deltaflag==1){pcell->marker.clone_index=1;}
		}
	}
	for(icell=1;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		if(nbneighb==5){
			deltaflag=1;
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){deltaflag=0;}
			}
			if(deltaflag==1){pcell->marker.clone_index=1;}
		}
	}
	for(icell=1;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		if(nbneighb==6){
			deltaflag=1;
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){deltaflag=0;}
			}
			if(deltaflag==1){pcell->marker.clone_index=1;}
		}
	}
}
int **rnd_parking(int *occupied, int *unoccupied, int *oclength, int *unoclength){
    Dual *pcell;
    int ocsize,unocsize,i,j,count=0;
    int index,cell,icell,delcount=0;
    int nbneighb;

    ocsize=*oclength;unocsize=*unoclength;
    int *tempoc1,*tempunoc1,*tempoc2,*tempunoc2;
    int *tmp;

    tempoc1=(int*) malloc(sizeof(int)*ocsize);
    if(tempoc1==NULL){locerror("rnd_parking","out of memory!");}
    tempunoc1=(int*) malloc(sizeof(int)*unocsize);
    if(tempunoc1==NULL){locerror("rnd_parking","out of memory!");}

    for(i=0;i<ocsize;i++){tempoc1[i]=occupied[i];}
    for(i=0;i<unocsize;i++){tempunoc1[i]=unoccupied[i];}

    //printf("we have ocsize = %d and unocsize = %d!\n",ocsize,unocsize);
    //sanity check!
    for(i=0;i<ocsize;i++){
        for(j=0;j<unocsize;j++){
            if(tempoc1[i]==tempunoc1[j]){locerror("rnd_parking","the vectors share a term!");}
        }
    }

    index=1+(int)(gsl_rng_uniform(rng)*(unocsize-1));
    if(TOROIDAL==1){index=(int)(gsl_rng_uniform(rng)*unocsize);}
    cell=tempunoc1[index];
    tempunoc1[index]=-1;delcount+=1;
    pcell=web_dual[cell];
    nbneighb=pcell->nb_vertices;
    for(i=0;i<nbneighb;i++){
        icell=indexcellinwebdual(pcell->celllist[i]);
        for(j=0;j<unocsize;j++){
            if(icell==tempunoc1[j]){tempunoc1[j]=-1;delcount+=1;}
        }
    }
    unocsize-=delcount;
    ocsize+=1;

    tempoc2=(int*) malloc(sizeof(int)*ocsize);
    if(tempoc2==NULL){locerror("rnd_parking","out of memory!");}
    tempunoc2=(int*) malloc(sizeof(int)*unocsize);
    if(tempunoc2==NULL){locerror("rnd_parking","out of memory!");}
    for(i=0;i<(ocsize-1);i++){tempoc2[i]=tempoc1[i];}
    tempoc2[ocsize-1]=cell;
    //printf("cell = %d, ocsize-1 = %d,tempoc2[ocsize-1]= %d",cell,ocsize-1,tempoc2[ocsize-1]);
    count=0;
    for(i=0;i<(unocsize+delcount);i++){
        if(tempunoc1[i]!=-1){tempunoc2[count]=tempunoc1[i];count++;}
    }
    if(count!=unocsize){
        printf("we have count = %d and unocsize = %d and delcount = %d\n",count,unocsize,delcount);
        locerror("rnd_parking()","something wrong with counting!");
    }

    int **result = malloc(2*sizeof(int*));
    result[0]=malloc(ocsize*sizeof(int));
    result[1]=malloc(unocsize*sizeof(int));
    for(i=0;i<ocsize;i++){result[0][i]=tempoc2[i];}
    for(i=0;i<unocsize;i++){result[1][i]=tempunoc2[i];}
    *oclength=ocsize;*unoclength=unocsize;

    free(tempoc1);free(tempoc2);
    free(tempunoc1);free(tempunoc2);

    return result;
}
void freeparking(int **park){
    int i;

    free(park[0]);
    free(park[1]);
    free(park);
}
void kagome_nn_modify(){ //places extra d-cells for cells completely surrounded by n-cells
    Dual *pcell,*pncell,*pnncell;
    int icell,jcell,kcell,deltaflag,nbneighb,count=0;
    int *occupied,*unoccupied,oclength,unoclength;
    int i,j,**result;

    //initialization in case there are cells already placed as deltas.
    oclength=0;unoclength=0;
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        nbneighb=pcell->nb_vertices;
        if(pcell->marker.clone_index==1){ //neighbors for delta clones
            oclength+=1;
            for(j=0;j<nbneighb;j++){
                pncell=pcell->celllist[j];
                if(pncell->marker.clone_index==0){pncell->marker.clone_index=2;}//set for convenience.  eliminated later
            }
        }
    }
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        if(pcell->marker.clone_index==0){unoclength+=1;}
    }
    occupied=(int*)malloc(oclength*sizeof(int));
    unoccupied=(int*)malloc(unoclength*sizeof(int));
    oclength=0;unoclength=0;
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        if(pcell->marker.clone_index==1){occupied[oclength]=i;oclength+=1;}
        if(pcell->marker.clone_index==0){unoccupied[unoclength]=i;unoclength+=1;}
    }
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        if(pcell->marker.clone_index==2){pcell->marker.clone_index=0;} //reset the index used for convenience
    }
    printf("placing delta cells...\n");

    while(unoclength>0){
        result=rnd_parking(occupied,unoccupied,&oclength,&unoclength);
        free(occupied);free(unoccupied);
        occupied=malloc(oclength*sizeof(int));unoccupied=malloc(unoclength*sizeof(int));
        for(i=0;i<oclength;i++){occupied[i]=result[0][i];}
        for(i=0;i<unoclength;i++){unoccupied[i]=result[1][i];}
        freeparking(result);
        /*
        printf("occupied elements: ");
        for(i=0;i<oclength;i++){printf("%d ",occupied[i]);}
        printf("\nunoccupied elements:");
        for(i=0;i<unoclength;i++){printf("%d ",unoccupied[i]);}
        printf("\n");fflush(stdout);
        */
    }
    printf("we have that %d cells are delta'd!\n",oclength);
    for(i=0;i<oclength;i++){
        pcell=web_dual[occupied[i]];
        pcell->marker.clone_index=1;
    }
}
void kagome_nnn_modify(){
   Dual *pcell,*pncell,*pnncell;
   int icell,jcell,kcell,deltaflag,nbneighb,count=0;

	printf("placing delta cells...\n");
    while(count<10000){
		deltaflag=1;icell=(int)(gsl_rng_uniform(rng)*nb_cell_tot);count++;
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		for(jcell=0;jcell<nbneighb;jcell++){
			pncell=pcell->celllist[jcell];
			if(pncell->marker.clone_index==1){deltaflag=0;}
            for(kcell=0;kcell<nbneighb;kcell++){
                pnncell=pncell->celllist[kcell];
                if(pnncell->marker.clone_index==1){deltaflag=0;}
            }
		}
		//if(deltaflag==1 && nbneighb>5){pcell->marker.clone_index=1;}
		if(deltaflag==1){pcell->marker.clone_index=1;}
	}
}
void kagome_modify2(){ //this is only for the specific toroid generated by alberto's code
    int i;
    Dual *pcell;

    printf("placing delta cells...\n");
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        if(i==1||i==7||i==13||i==19||i==26||i==32||i==38||i==44||i==4||i==10||i==16||i==22||i==29||i==35||i==41||i==47){pcell->marker.clone_index=1;}
    }
}
void kagome_dblmodify(){
	Dual *pcell,*pncell;
	int icell,jcell,deltaflag,nbneighb,count,numdelta=0,numnotch=0;
	int numchange;

	printf("placing secondary delta cells...\n");
	for(icell=0;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		deltaflag=0;count=0;
		if(pcell->marker.clone_index==1){numdelta++;}
		if(pcell->marker.clone_index==0){numnotch++;}
		if(nbneighb>4 && pcell->marker.clone_index==0){
			deltaflag=1;
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){count++;}
			}
			if(count>1){deltaflag=0;}
		}
		if(deltaflag==1){pcell->marker.clone_index=1;numdelta++;numnotch--;}
	}
	for(icell=0;icell<nb_cell_tot;icell++){ //this is to eliminate the rare double event
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		count=0;
		if(pcell->marker.clone_index==1){
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){count++;}
			}
		}
		if(count>1){pcell->marker.clone_index=0;numdelta--;numnotch++;}
	}
	printf("we have %d delta cells and %d notch cells!\n",numdelta,numnotch+1);fflush(stdout);
	if(numdelta>numnotch/2-nb_cell_tot/80){numchange=(int)ceil(numdelta-numnotch/2+nb_cell_tot/80);}
	else{numchange=0;}
	printf("we need to change %d cells!\n",numchange);fflush(stdout);
	icell=1; //initialize
	while(numchange>0){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;count=0;
		for(jcell=0;jcell<nbneighb;jcell++){
			pncell=pcell->celllist[jcell];
			if(pncell->marker.clone_index==1){count++;}
		}
		if(count>0 && pcell->marker.clone_index==1){pcell->marker.clone_index=0;numchange--;printf("changing cell %d...",icell);}
		icell++;
		if(icell>nb_cell_tot-1){break;}
	}
	printf("done!\n");fflush(stdout);
}
void kagome_reg(){
    Dual *pcell,*pncell;
    int i,j,nbneighb,delta;

    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        delta=0;
        if(pcell->marker.clone_index==0){
            delta=1;nbneighb=pcell->nb_vertices;
            for(j=0;j<nbneighb;j++){
                pncell=pcell->celllist[j];
                if(pncell->marker.clone_index==1){delta=0;}
            }
        }
        if(delta==1){pcell->marker.clone_index=1;}
    }
}

   void
init_lattice2(int *array, int num)
{
   int icell,j,k;
   Dual *pcell;
   char filename1[100],filename2[100],fnum[4];


   printf("before coord\n");
   //coord_gen();//comment this out when using infoinput()
   printf("before shrink\n");
   sprintf(fnum,"%.4d",array[num]);
   sprintf(filename1,"info/%s.vertinfo",fnum);
   sprintf(filename2,"info/%s.cellinfo",fnum);
   infoinput(filename1,filename2,0);
   //cut_cell_find();
   //shrink_datafirst();//comment this out to use infoinput()
   //coord_mod();//to give the initial vertices an offset
   printf("after shrink\n");


   for (icell=0 ; icell<nb_cell_tot ; icell++){
      pcell=web_dual[icell];
      pcell->marker.size=0;
      pcell->marker.div_rate=0;
      pcell->marker.clone_index=0;
      pcell->area_soll=1.0;
      pcell->area_soll0=pcell->area_soll;
      pcell->marker.tension_index=1.0;
   }

   forcevector=(double**)malloc(3*sizeof(double*));
   for(k=0;k<3;k++){
      forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
   }
   //svglattice("svgs/0000.svg",1);
   //latticeprint("dats/0000.dat",0);
   //cut_cell_multistep();
   /*for(icell=0;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      if(icell==20 || icell==23 || icell==27 || icell==24 || icell==30){pcell->marker.tension_index+=add_factor;}
      if(icell==0){pcell->area_soll=78;pcell->area_soll0=pcell->area_soll;}
      optimizet1();
   }*/
   itime = 0;
   it1process=0;
}

   void
import_lattice(char * filename)
{
   int icell,j,k;
   Dual *pcell;
   char filename1[100],filename2[100];


   sprintf(filename1,"input/%s.vertinfo",filename);
   sprintf(filename2,"input/%s.cellinfo",filename);
   infoinput(filename1,filename2,1);

   for (icell=0 ; icell<nb_cell_tot ; icell++){
      pcell=web_dual[icell];
      pcell->marker.size=0;
      pcell->marker.div_rate=0;
      pcell->area_soll=pcell->area_soll;
      pcell->area_soll0=pcell->area_soll;
      pcell->marker.tension_index=pcell->area_soll;
   }

   forcevector=(double**)malloc(3*sizeof(double*));
   for(k=0;k<3;k++){
      forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
   }
   itime = 0;
   it1process=0;
}

void kagomet1_dd(){ //transforms delta cells away from each other
	int i,j,k,l,nbneighb,exec_t1;
	Dual *pcell,*pncell;
	Node *pvert,*pnvert,*povert;
	int debugcount=0,start;
	char filename[50];

    start=1;
    if(TOROIDAL==1){start=0;}
	//this first part does cells for neighboring delta'd cells
	for(i=start;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		if(pcell->marker.clone_index==1){
			//printf("we found a cell to transform! it is cell %d\n",i);fflush(stdout);
			nbneighb=pcell->nb_vertices;
			//here we check to see if any delta'd cells are neighbors.  if so, then transform them apart. assumes at most 1 neighboring delta.
			exec_t1=0;
			for(j=0;j<nbneighb;j++){
				pncell=pcell->celllist[j];
                //printf("we are in cell %d ",indexcellinwebdual(pcell));fflush(stdout);
                //printf("and neighbor %d ",j);fflush(stdout);
                //printf("which has global index %d!\n",indexcellinwebdual(pncell));fflush(stdout);
				if(pncell->marker.clone_index==1){
					//printf("okay we found two deltas next to each other. cell %d and neighbor %d...",i,j);fflush(stdout);
					//this part just finds the index of the other vertex in the pneighb list on vertex structures
					exec_t1=1;
					pvert=pcell->vertexlist[j];pnvert=pcell->vertexlist[(j+1)%nbneighb];//transform vertices on the boundary between delta'd cells
					for(k=0;k<3;k++){
						povert=pvert->pneighb[k];
						if(povert==pnvert){l=k;}
					}
				}
			}
			if(exec_t1==1){t1transform(pvert,l);/*printf("transformed!\n");fflush(stdout);*/} //transform to separate the delta cells
		}
	}
	updatecellbondneighbs();
}
void kagomet1_nd(){
	int i,j,k,l,nbneighb,exec_t1;
	Dual *pcell,*pncell;
	Node *pvert,*pnvert,*povert;
	int debugcount=0,start;
	char filename[50];

    start=1;
    if(TOROIDAL==1){start=0;}
	for(i=start;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		if(pcell->marker.clone_index==1){
			//printf("we found a cell to transform! cell %d\n",i);fflush(stdout);
			//here we try to stick as many notch to the deltas as we can thru t1transforms
			nbneighb=pcell->nb_vertices;
			for(j=0;j<nbneighb;j++){
				//printf("looking at vertex %d in cell %d...",j,i);fflush(stdout);
				pvert=pcell->vertexlist[j];
				for(k=0;k<3;k++){
					if((pvert->pneighb[k]!=pcell->vertexlist[(j-1)%nbneighb]) && (pvert->pneighb[k]!=pcell->vertexlist[(j+1)%nbneighb])){pnvert=pvert->pneighb[k];l=k;}
				}
				exec_t1=1;
				for(k=0;k<3;k++){//this determines whether or not any neighbors are delta'd. if so, then don't execute t1
					pncell=pnvert->pncell[k];
					if(pncell->marker.clone_index==1){
						//printf("awww, no transform :(\n");fflush(stdout);
						exec_t1=0;
					}
				}
				if(exec_t1==1){
					//printf("performing t1!\n");fflush(stdout);
					t1transform(pvert,l);
					optimize();
					sprintf(filename,"debugsvg/%.4i.svg",debugcount);
					//svglattice(filename,1);//uncomment to see each step
					debugcount++;
				}
			}
		}
	}
	updatecellbondneighbs();

    Bond *pbond;
    int m;

    for(i=0;i<nb_bond_tot;i++){
        pbond=web_bond[i];m=0;l=0;k=-1;
        pvert=pbond->pnvert[0];pnvert=pbond->pnvert[1];
        for(j=0;j<3;j++){
            if(pvert->pneighb[j]==pnvert){k=j;}
        }
        if((pbond->pncell[0])->marker.clone_index==0){l+=1;}
        if((pbond->pncell[2])->marker.clone_index==0){l+=1;}
        if((pbond->pncell[1])->marker.clone_index==0){m+=1;}
        if((pbond->pncell[3])->marker.clone_index==0){m+=1;}
        if(l==2 && m==1){t1transform(pvert,k);}
    }
    updatecellbondneighbs();
}

void pulse_event(double step_frac,double n_str){

	int i,j;
	Dual *pcell;

	for(i=1;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		pcell->area_soll=pcell->area_soll-pcell->area_soll0*step_frac+2*n_str*(gsl_rng_uniform(rng)-0.5);
	}
}

int *division_select_once(int divlist[],int divlistlen,int *pnum){ //selects one number from array and removes it

	int i,idx,len,*templist;

	len=divlistlen;
	idx=floor(gsl_rng_uniform(rng)*len);//select one cell randomly
	*pnum=divlist[idx];

	//temporarily store the remaining elements
	templist=(int*)malloc((len-1)*sizeof(int));
	for(i=0;i<idx;i++){templist[i]=divlist[i];}
	for(i=idx;i<len-1;i++){templist[i]=divlist[i+1];}

	return templist;
}

void torus_division_set_parameters(){

	int i,j;
	Dual *pcell;

	for(i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];

		pcell->marker.tension_index=3.0;
		pcell->area_soll=LX*LY/nb_cell_tot;
		pcell->area_soll0=pcell->area_soll;
	}
}

Dual *findCenterCell(){
	int i,idx;
	Dual *pcell;
	double **pos;
	vector3 cent,center;

	center.x=0.;center.y=0.;center.z=0.;
	pos=(double**)calloc(3,sizeof(double*));
	for(i=0;i<3;i++){pos[i]=(double*)calloc(nb_cell_tot,sizeof(double));}
	for(i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		cent=pcell->centroid;
		pos[0][i]=cent.x;pos[1][i]=cent.y;
		if(NUMDIM==3) pos[2][i]=cent.z;
		center.x+=(pos[0][i]/((double)nb_cell_tot));
		center.y+=(pos[1][i]/((double)nb_cell_tot));
		center.z+=(pos[2][i]/((double)nb_cell_tot));
	}
	double dx,dy,dz,d;
	double mindist=100.0;
	for(i=0;i<nb_cell_tot;i++){
		dx=center.x-pos[0][i];dy=center.y-pos[1][i];
		if(NUMDIM==3){dz=center.z-pos[2][i];}
		else if(NUMDIM==2){dz=0.;}
		d=sqrt(dx*dx+dy*dy+dz*dz);
		if(d<mindist){mindist=d;idx=i;}
	}
	printf("we found cell %d to be closest to the center!\n",idx);

	for(i=0;i<NUMDIM;i++){free(pos[i]);}
	free(pos);
	return web_dual[idx];
}


gsl_matrix *computeDynMat(){ //compute the dynamical matrix

	int i,j,k;
	double A,A0,L,P,P0,E=0.,kP,kA;
	double dE_dA,dE_dL,dE_dP;
	double n_area,n_perm,n_length;
	Dual *pcell;

	////the following code block was to check energy. outdated.////
	/*
	//cell contributions
	for(i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		A=pcell->area;P=pcell->perimeter;A0=pcell->area_soll;P0=pcell->perim_soll;
		E+=kA*(A-A0)*(A-A0)+kP*(P-P0)*(P-P0);
	}
	//bond contributions
	for(i=0;i<nb_bond_tot;i++){
		//none for now
		//E+=0.;
		break;
	}

	//make sure the model matches relaxation performed by opti.c
	double Etest=0.;
	for(i=0;i<nb_cell_tot;i++){
		Etest+=kP*rho_fct(pcell)+kA*area_fct(pcell)+extern_pot_fct(pcell);
	}
	if(fabs(Etest-E)/Etest > 1e-10) locerror("computeDynMat()","Dynamical Matrix energy form does not agree with optimization scheme!\n");
	*/

	//Consider a d=2 model for now
	int dim=NUMDIM;

	//allocate the matrix
	gsl_matrix *dynmat;
	double *dPdr1,*dPdr2,**d2Pdr1dr2,*temp; //perimeter terms
	double *dAdr1,*dAdr2,**d2Adr1dr2;
	int idx1,idx2;
	Node *pv1,*pv2;

	//kP is scaled out in the computation of the dynamical matrix for simplicity

	//dynamical matrix memory allocation
	dynmat=gsl_matrix_calloc(dim*nb_vertex_tot,dim*nb_vertex_tot);
	//Perimeter term contributions: memory allocation
	dPdr1=(double*)calloc(dim,sizeof(double));dPdr2=(double*)calloc(dim,sizeof(double));
	d2Pdr1dr2=(double**)calloc(dim,sizeof(double*));
	for(i=0;i<dim;i++){d2Pdr1dr2[i]=(double*)calloc(dim,sizeof(double));}
	//Area term contributions: memory allocation
	dAdr1=(double*)calloc(dim,sizeof(double));dAdr2=(double*)calloc(dim,sizeof(double));
	d2Adr1dr2=(double**)calloc(dim,sizeof(double*));
	for(i=0;i<dim;i++){d2Adr1dr2[i]=(double*)calloc(dim,sizeof(double));}

	for(i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		//compute mode contributions from this cell vertex by vertex
		for(j=0;j<pcell->nb_vertices;j++){
			pv1=pcell->vertexlist[j];idx1=pv1->idx;
			dP_dr(pcell,j,dPdr1);
			dA_dr(pcell,j,dAdr1);
			for(k=0;k<pcell->nb_vertices;k++){
				pv2=pcell->vertexlist[k];idx2=pv2->idx;
				dP_dr(pcell,k,dPdr2);
				dA_dr(pcell,k,dAdr2);

				//The following lines add on contributions from perimeter terms
				//quadratic term for perimeter: kP*dL_j*dL_k
				outerProd(dPdr1,dPdr2,d2Pdr1dr2);
				//This term has a prefactor of second derivative d2E_dP2
				addIntoDynMat(dynmat,d2Pdr1dr2,idx1,idx2,d2E_dP2(pcell));
				//check whether or not pv1 and pv2 share a bond or are same vertex
				//if it does, then add second order terms which have a prefactor 
				//involving the pre-stress:dE_dP (rho_fctpp())
				if(d2P_drdr(pcell,j,k,d2Pdr1dr2)==1){addIntoDynMat(dynmat,d2Pdr1dr2,idx1,idx2,rho_fctpp(pcell));}


				//The following lines add on contributions from area terms
				//quadratic term for area: dA*dA_j*dA_k
				outerProd(dAdr1,dAdr2,d2Adr1dr2);
				//this term has a prefactor of second derivative d2E_dA2
				addIntoDynMat(dynmat,d2Adr1dr2,idx1,idx2,d2E_dA2(pcell));
				//check whether or not pv1 and pv2 share a bond (cannot be same vertex)
				//if they do, add second order (levi-civita) term. This term has a
				//prefactor involving the pre-stress: dE_dA (area_fctp())
				if(d2A_drdr(pcell,j,k,d2Adr1dr2)==1){addIntoDynMat(dynmat,d2Adr1dr2,idx1,idx2,area_fctp(pcell));}
			}
		}


	}

	free(dPdr1);free(dPdr2);free(dAdr1);free(dAdr2);
	for(i=0;i<NUMDIM;i++){free(d2Pdr1dr2[i]);free(d2Adr1dr2[i]);}
	free(d2Pdr1dr2);free(d2Adr1dr2);

	return dynmat;
}

int outerProd(double *vec1,double *vec2,double **outmat){
	int i,j;
	for(i=0;i<NUMDIM;i++){for(j=0;j<NUMDIM;j++){outmat[i][j]=vec1[i]*vec2[j];}}
	return 1;
}

int addIntoDynMat(gsl_matrix *dynmat,double **inmat,int idx1,int idx2,double prefactor){
	int i,j;
	double tempval;
	for(i=0;i<NUMDIM;i++){for(j=0;j<NUMDIM;j++){
		tempval=gsl_matrix_get(dynmat,NUMDIM*idx1+i,NUMDIM*idx2+j);
		gsl_matrix_set(dynmat,NUMDIM*idx1+i,NUMDIM*idx2+j,tempval+inmat[i][j]);
	}}
	//scale everything up by input prefactor
	gsl_matrix_scale(dynmat,prefactor);
	return 1;
}

int issameMatrices(const gsl_matrix *mat1,const gsl_matrix *mat2,int dim1,int dim2){
	int i,j,num;
	gsl_matrix *outmat;
	double temp1,temp2;
	int count=0;
	double tol=1e-10; //error tolerance between check values

	//variables to track entries that are different
	double *vals;
	int *indices;

	//not even the same size
	if(mat1->size1!=mat2->size1 || mat1->size2!=mat2->size2){return -1;}
	
	dim1=mat1->size1;dim2=mat1->size2;

	outmat=gsl_matrix_calloc(dim1,dim2);
	gsl_matrix_memcpy(outmat,mat1);
	gsl_matrix_sub(outmat,mat2); //subtracts mat2 from outmat (copied form mat1)
	for(i=0;i<dim1;i++){
		for(j=0;j<dim2;j++){
			temp1=fabs(gsl_matrix_get(outmat,i,j)/gsl_matrix_get(mat1,i,j));
			temp2=fabs(gsl_matrix_get(outmat,i,j)/gsl_matrix_get(mat2,i,j));
			//find largest difference
			if(temp2>temp1){gsl_matrix_set(outmat,i,j,temp2);}
			else{gsl_matrix_set(outmat,i,j,temp1);}
			//compare to tolerance
			if(gsl_matrix_get(outmat,i,j)>tol){count++;}
		}
	}
	if(count==0){ //return success! else continue for diagnostics
		gsl_matrix_free(outmat);
		return 1;
	}

	num=count;count=0;
	vals=(double*)calloc(num,sizeof(double));
	indices=(int*)calloc(2*num,sizeof(int));
	//create data arrays of values that exceed tolerance
	for(i=0;i<dim1;i++){for(j=0;j<dim2;j++){
		temp1=gsl_matrix_get(outmat,i,j);
		if(temp1>tol){vals[count]=temp1;indices[2*count]=i;indices[2*count+1]=j;count++;}
	}}
	printf("The following entries are unequal beyond tolerance:\n");
	for(i=0;i<num;i++){printf("Difference between matrices is %lf at (%d,%d)\n",vals[i],indices[2*i],indices[2*i+1]);}

	gsl_matrix_free(outmat);
	free(vals);free(indices);
	return 0;
}

int addVal_gsl_matrix(gsl_matrix *m,size_t i,size_t j,double val){
	gsl_matrix_set(m,i,j,gsl_matrix_get(m,i,j)+val);
	return 1;
}

gsl_matrix *computePerimCompatMat(){
	int i,j,k;
	Dual *pcell;
	Node *pv1,*pv2;
	gsl_matrix *C;
	double *tempvec;
	double fac;

	//Just like in the regular dynamical matrix, we normalize by
	//kP, the perimeter pre-factor
	C=gsl_matrix_calloc(NUMDIM*nb_vertex_tot,nb_cell_tot+nb_bond_tot);
	//split it into the following types of "extension": total
	//parallel perimeter and bond-wise perpendicular

	//data vector to take values temporarily
	tempvec=(double*)calloc(NUMDIM,sizeof(double));
	for(i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		for(j=0;j<pcell->nb_vertices;j++){
			pv1=pcell->vertexlist[j];pv2=pcell->vertexlist[(j+1)%pcell->nb_vertices];
			parallelProj(pv1,pv2,tempvec);
			addVal_gsl_matrix(C,i,NUMDIM*(pv1->idx),tempvec[0]);
			addVal_gsl_matrix(C,i,NUMDIM*(pv1->idx)+1,tempvec[1]);
			addVal_gsl_matrix(C,i,NUMDIM*(pv2->idx),-1.*tempvec[0]);
			addVal_gsl_matrix(C,i,NUMDIM*(pv2->idx)+1,-1.*tempvec[1]);
			if(NUMDIM==3){
				addVal_gsl_matrix(C,i,NUMDIM*(pv1->idx)+2,tempvec[2]);
				addVal_gsl_matrix(C,i,NUMDIM*(pv2->idx)+2,-1.*tempvec[2]);
			}
		}
	}

	Dual *pcell1,*pcell2;

	for(i=0;i<nb_bond_tot;i++){
		j=i+nb_cell_tot; //re-index for each bond-based extension
		pcell1=web_bond[i]->pncell[0];pcell2=web_bond[i]->pncell[2];
		pv1=web_bond[i]->pnvert[0];pv2=web_bond[i]->pnvert[1];
		//This factor is T / R -- tension divided by the bond length
		fac=((pcell1->perimeter-pcell1->perim_soll)+(pcell2->perimeter-pcell2->perim_soll))/dist(pv1,pv2);
		d2RotatePerpProj(pv1,pv2,tempvec); //do this in d=2
		//d=2 is special, since there is only one perpendicular direction we can just
		//effectively represent it by e^a = \epsilon^{ab} R^b/R and therefore easily
		//find the magnitude of \delta R^\perp. Using the usual transverse projector
		//\delta^{ab}-R^a R^b / R^2 is like defining the sin in terms of the cosine:
		//perfectly reasonable but no longer linear in the argument (cos = sqrt(1-sin^2))
		//which (at least potentially) creates issues when constructing the compatibility matrix

		//we do not incorporate a d=3 version as a result. The counting is going to be very
		//different anyways, and perhaps it may not be relevant for these systems we hope
		//to study that are near Maxwell. However, d=3 for the dynamical
		//matrix may still be interesting, so I will consider coding the rest of that up.
		addVal_gsl_matrix(C,j,NUMDIM*(pv1->idx),sqrt(fac)*tempvec[0]);
		addVal_gsl_matrix(C,j,NUMDIM*(pv1->idx)+1,sqrt(fac)*tempvec[1]);
		addVal_gsl_matrix(C,j,NUMDIM*(pv2->idx),-1.*sqrt(fac)*tempvec[0]);
		addVal_gsl_matrix(C,j,NUMDIM*(pv2->idx)+1,-1.*sqrt(fac)*tempvec[1]);
		if(NUMDIM==3){
			//have an extra transverse direction
		}
	}

	free(tempvec);
}
