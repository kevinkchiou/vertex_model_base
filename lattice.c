/* lattice3-1-04
 * Evolution of a lattice with cell divisions
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>

#include "const.h" /* Main file where the parameters are defined */

#include "definelat.h"
#include "generatelat.h"
#include "save.h"
#include "read.h"
#include "lattice.h" /* Declaration of the structures of the lattice and declaration of the functions defined in this file */
#include "locerror.h" /* Error function declaration */
#include "measurements.h" /* Declaration of the functions used for statistics */
#include "pngwrite.h"
#include "opti.h"
#include "compalone.h"
#include "in_development.h"



/* Gives the indice of a vertix from the cartesian coordinates on the first lattice
*/
   int
vertex_from_coord(int row, int column)
{
   if (row>2*LENGTH_FIRST_LATTICE-1 || column>LENGTH_FIRST_LATTICE){
      locerror("vertex_from_coord", "coordinates went over the limits");
   }
   return LENGTH_FIRST_LATTICE*(row-1)+column;
}


/* Initializes the first regular lattice
*/
   void
init_lattice()
{
	int icell,j,k;
	Dual *pcell;
	char filename1[100],filename2[100];


	printf("before coord\n");
	coord_gen();//comment this out when using infoinput()
	printf("before shrink\n");
	sprintf(filename1,"info/0100.vertinfo");
	sprintf(filename2,"info/0100.cellinfo");
	//infoinput(filename1,filename2);
	//cut_cell_find();
	shrink_datafirst();//comment this out to use infoinput()
	//coord_mod();//to give the initial vertices an offset
	printf("after shrink\n");

	for(k=0;k<nb_vertex_tot;k++){web[k]->idx=k;}

	for (icell=0 ; icell<nb_cell_tot ; icell++){
		pcell=web_dual[icell];
		pcell->marker.size=0;
		pcell->marker.div_rate=0;
		pcell->marker.clone_index=0;
		pcell->area_soll=1;
		pcell->area_soll0=pcell->area_soll;
		pcell->sigma=1.0;
	}

	forcevector=(double**)malloc(3*sizeof(double*));
	for(k=0;k<3;k++){
		forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
	}
	svglattice("svgs/0000.svg",1);
	latticeprint("dats/0000.dat",0);
	//cut_cell_multistep();
	createcellbondneighbs();optimizet1();
}

void init_lattice_torus()
{
	int n,i,k,l,m,kk;
	double sq3=sqrt(3.),xc,yc,cox,coy,tol=0.001,LXLOC,LYLOC;
	Dual *pcell;

	//coord_gen_torus2();
	coord_gen_torus();
	for (i=0 ;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		pcell->marker.size=0;
		pcell->marker.div_rate=0;
		pcell->marker.clone_index=0;

		//put the tissue under tensile stress
		//pcell->area_soll0=(LX*LY)/nb_cell_tot; //original A0. should not change
		//pcell->area_soll=pcell->area_soll0; //A0 in energy. can change
		//pcell->marker.tension_index=pcell->area_soll;

		//put the tissue under compressive stress
		pcell->area_soll0=(LX*LY*2.)/nb_cell_tot;
		pcell->area_soll=pcell->area_soll0;
		pcell->marker.tension_index=1.0;

		pcell->marker.border=0;
	}
	forcevector=(double**)malloc(3*sizeof(double*));
	for(k=0;k<3;k++){
		forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
	}
	createcellbondneighbs();
	optimizet1();

}

/* Free the memory used to store the lattice */
	void
kill_lattice()
{
	int i;

	cleanup(); //clean up  data structures with redundant geometric information (bonds, cell neighbors)
	for (i=0; i<nb_cell_tot; i++){
		freecell(web_dual[i]);web_dual[i]=NULL;
	}

	for (i=0; i<nb_vertex_tot; i++){
		free(web[i]);web[i]=NULL;
	}

	free(web);web=NULL;
	free(web_dual);web_dual=NULL;
}

/* If a pointer to a cell or a vertex is void, the array is shrinked so that at the end the number of pointers is equal to the number of cells of vertices.
*/
void shrink_datafirst()
{
	int ivertex, idefvert, icell, idefcell;
	Node **webtemp;
	Dual **web_dualtemp;
	/* Node list shrinking */

	nb_vertex_tot=0;

	for (ivertex=0 ; ivertex < last_first_lattice ; ivertex++){
		if (web[ivertex]){
			nb_vertex_tot++;
		}
	}
	webtemp=createnodevector(nb_vertex_tot);

	//printf("Total number of vertices : %i\n", nb_vertex_tot);

	idefvert=0;

	for (ivertex=0 ; ivertex < last_first_lattice ; ivertex++){
		if (web[ivertex]){
			webtemp[idefvert]=web[ivertex];
			idefvert++;
		}
	}
	free(web);
	web=webtemp;

	/* Cell list shrinking */    
	nb_cell_tot=0;

	for (icell=0 ; icell<last_first_lattice ; icell++){
		if (web_dual[icell]){
			nb_cell_tot++;
		}
	}
	web_dualtemp=createdualvector(nb_cell_tot);

	//printf("Total number of cells : %i\n", nb_cell_tot);fflush(stdout);

	idefcell=0;

	for (icell=0 ; icell<last_first_lattice ; icell++){
		if (web_dual[icell]){
			web_dualtemp[idefcell]=web_dual[icell];
			idefcell++;
		}
	}
	free(web_dual);
	web_dual=web_dualtemp;
}

/* Shrink web[] and web_dual[] when a cell is being killed
*/
void shrink_data()
{
	int ivertex, icell;
	int nb_vertex_def=0;
	int nb_cell_def=0;
	int idefvert=0;
	int idefcell=0;
	Node **webtemp;
	Dual **web_dualtemp;

	/* Node list shrinking */
	for (ivertex=0 ; ivertex < nb_vertex_tot ; ivertex++){
		if (web[ivertex]){
			nb_vertex_def++;
		}
	}
	webtemp=createnodevector(nb_vertex_def);

	for (ivertex=0 ; ivertex < nb_vertex_tot ; ivertex++){
		if (web[ivertex]){
			webtemp[idefvert]=web[ivertex];
			idefvert++;
		}
	}
	nb_vertex_tot=nb_vertex_def;
	//printf("Total number of vertices : %d\n",nb_vertex_tot);

	free(web);
	web=webtemp;

	/* Cell list shrinking */    
	for (icell=0 ; icell<nb_cell_tot ; icell++){
		if (web_dual[icell]){
			nb_cell_def++;
		}
	}
	web_dualtemp=createdualvector(nb_cell_def);

	for (icell=0 ; icell<nb_cell_tot ; icell++){
		if (web_dual[icell]){
			web_dualtemp[idefcell]=web_dual[icell];
			idefcell++;
		}
	}
	nb_cell_tot=nb_cell_def;
	//printf("Total number of cells : %d\n",nb_cell_tot);fflush(stdout);

	free(web_dual);
	web_dual=web_dualtemp;
}

/* Computes the exterior cell composed of the vertices on the boundary
*/
void cell0gen(Node *pfirstvertex)
{
	int icell, vertcount;
	Node *pvertex;
	Dual *pcell0;

	pvertex=NULL;
	vertcount=1;
	pcell0=web_dual[0];

	icell=indexcell(pfirstvertex, pcell0);
	pvertex=pfirstvertex -> pneighb[icell];
	while(pvertex!=pfirstvertex){
		vertcount++;
		icell=indexcell(pvertex, pcell0);
		pvertex=pvertex -> pneighb[icell];
	}

	free(pcell0 -> vertexlist);
	pcell0 -> vertexlist = createnodevector(vertcount);
	pcell0 -> nb_vertices = vertcount;

	pvertex=NULL;
	vertcount=0;
	pcell0 -> vertexlist[0]=pfirstvertex;
	icell=indexcell(pfirstvertex, pcell0);
	pvertex=pfirstvertex -> pneighb[icell];
	while(pvertex!=pfirstvertex){
		vertcount++;
		icell=indexcell(pvertex, web_dual[0]);
		pvertex=pvertex -> pneighb[icell];
	}
}




/* indexcell return the index of cell for a given vertex and cell */ 
	int
indexcell(Node *pvertex, Dual *pcell)
{
	int i;

	for (i=0; i<3 ; i++ ){
		if ((pvertex -> pncell[i])==pcell) return i;
	}
	locerror("indexcell", "the cell does not belong to the vertex");
	return -1;
}

int indexvertex(Node *pvertex, Dual *pcell)
{
	int nbvertices, n;

	nbvertices=pcell->nb_vertices;

	for (n=0; n < nbvertices ; n++){
		if ((pcell->vertexlist[n])==pvertex){
			return n;
		}
	}
	locerror("indexvertex", "The vertex does not belong to the cell");
	return -1;
}

	int
indexvertexinweb(Node *pvertex)
{
	int n;

	for (n=0; n < nb_vertex_tot ; n++){
		if (web[n]==pvertex){
			return n;
		}
	}
	locerror("indexvertexinweb", "The vertex does not belong to the web!!");
	return 0;
}

	int
indexcellinwebdual(Dual *pcell)
{
	int n;

	for (n=0; n < nb_cell_tot ; n++){
		if (web_dual[n]==pcell){
			return n;
		}
	}
	locerror("indexcellinwebdual", "The cell does not belong to the web_dual!!");
	return 0;
}

/* Insert a new vertex in vertexlist */
	void
insvertex(Dual *pcell, Node *pnewvert, Node *pnextvert) /* The new vertex is at the first position */
{
	int nbvertices, n, ivert;
	Node **pvertexlisttemp;

	if ((ivert=indexvertex(pnextvert, pcell)) < 0){
		printf("In insvertex ivert < 0\n");
		exit(1);
	};

	nbvertices=pcell->nb_vertices;
	pvertexlisttemp=createnodevector(nbvertices+1);
	pcell->nb_vertices =nbvertices+1;

	pvertexlisttemp[0]=pnewvert;
	for (n=0 ; n < nbvertices ; n++){
		pvertexlisttemp[n+1]=pcell->vertexlist[(ivert+n)%nbvertices];
	}

	free(pcell->vertexlist);
	pcell->vertexlist = pvertexlisttemp;
}

int indexbondinwebbond(Bond *pbond){
	int i;

	for(i=0;i<nb_bond_tot;i++){
		if(web_bond[i]==pbond){return i;}
	}
	locerror("indexbondinwebbond","The bond does not exist");
	return 0;
}

	void
delvertex(Dual *pcell, Node *pdelvert)
{
	int ivert, n, nbvertices;
	Node **pvertexlisttemp;

	if ((ivert=indexvertex(pdelvert, pcell)) < 0){
		printf("In delvertex ivert < 0\n");
		exit(1);
	};


	nbvertices=pcell->nb_vertices;
	pvertexlisttemp=createnodevector(nbvertices-1);
	pcell->nb_vertices =nbvertices-1;

	for (n=0 ; n < nbvertices-1 ; n++){
		pvertexlisttemp[n]=pcell->vertexlist[(ivert+n+1)%nbvertices];
	}

	free(pcell->vertexlist);
	pcell->vertexlist = pvertexlisttemp;
}


	int
kill_cell_p(Dual *pcell)
{
	int i, j, ix;
	Node *ptmpvert,*ptttvert;
	Node *cv[2];
	Dual *pcelltmp,*pcell0, *pcellp;
	double pres, presdif;


	pcell0 = pcellp = web_dual[0];

	// find neighboring cell with highest/lowest pressure
	pres = pressure(pcell);
	printf("%f\n",pres);
	presdif = 0.;
	for (i=0;i<pcell->nb_vertices;i++) {
		ptmpvert = pcell->vertexlist[i];
		for (j=0;j<3;j++){
			pcelltmp = ptmpvert->pncell[j];
			if (pcelltmp != pcell && pcelltmp != pcell0){
				printf("%i %i %f\n",i,j,pressure(pcelltmp));
				if (fabs(pressure(pcelltmp)-pres) > presdif) {
					pcellp  = pcelltmp;
					presdif = fabs(pressure(pcelltmp)-pres);
				}
			}
		}
	}
	printf("%f %f\n",pres, pressure(pcellp));

	if (pcell->nb_vertices == 3 && pcellp->nb_vertices == 3) {
		printf("pcell and pcellp have both 3 vertices\n");
		exit(1);
	}

	// find verticies commen two pcell and pcellp
	printf("pcell  has %i vertices\n",pcell->nb_vertices);
	printf("pcellp has %i vertices\n",pcellp->nb_vertices);
	ix = 0;
	for (i=0;i<pcell->nb_vertices;i++) {
		ptmpvert = pcell->vertexlist[i];
		for (j=0;j<pcellp->nb_vertices;j++) {
			ptttvert = pcellp->vertexlist[j];
			if (ptttvert == ptmpvert) {
				if (ix >1) {
					printf("pcell and pcellc hav more than 2 vertices in commen\n");
					exit(1);
				}
				cv[ix] = ptttvert;
				ix++;
			}
		}
	}

	// check neighboring cells
	for (j=0;j<3;j++){
		ptmpvert = cv[0];
		pcelltmp = ptmpvert->pncell[j];
		if (pcelltmp != pcell && pcelltmp != pcell0 && pcelltmp->nb_vertices == 3){
			printf("cv[1] neighb cell has both 3 vertices\n");
			exit(1);
		}
		ptmpvert = cv[1];
		pcelltmp = ptmpvert->pncell[j];
		if (pcelltmp != pcell && pcelltmp != pcell0 && pcelltmp->nb_vertices == 3){
			printf("cv[2] neighb cell has both 3 vertices\n");
			exit(1);
		}
	}

	for (j=0;j<3;j++) {
		ptmpvert = cv[0]->pneighb[j];
	}



	exit(1);
}


int find_smallc(Dual *pcell)
{
	int nvert, index=0;
	double areamin=100000;

	for (nvert=0;nvert<pcell->nb_vertices;nvert++){
		if (pcell->area<areamin){
			areamin=pcell->area;
			index=nvert;
		}
	}
	return index;
}

	int
kill_cell_area(Dual *pcell)
{
	int index;
	vector3 center;

	index=find_smallc(pcell);
	center=centroid(pcell);
	kill_cell(pcell,index);
	return 0;
}

/* Kill a cell */
	int
kill_cell(Dual *pcell, int ivertex)
{
	int i, icell, icell1, icell2, icella, icellb, icellc, icelld, ivertexbig;
	int nbvertices1, nbvertices2, clone;
	//char filename[256];
	Node *pvert;
	Node *pvertex1, *pvertex2, *pvertexa, *pvertexb, *pvertexc, *pvertexd;
	Node **vertexlisttemp;
	Dual *pcellbig,*pcelltmp,*pcell0;
	double pres, presdif;


	pcell0 = web_dual[0];
	clone = 0;
	if (pcell->marker.div_rate == 1) {
		clone = 1;
	}

	// find neighboring cell with highest/lowest pressure
	pres = pressure(pcell);
	//printf("%f\n",pres);
	presdif = 0.;
	/*for (i=0;i<pcell->nb_vertices;i++) {
	  ptmpvert = pcell->vertexlist[i];
	  j=indexcell(ptmpvert, pcell);
	  pcelltmp=ptmpvert->pncell[(j+2)%3];

	  if (pcelltmp != pcell && pcelltmp != pcell0){
	//printf("%i %i %f\n",i,j,pressure(pcelltmp));
	if (fabs(pressure(pcelltmp)-pres) > presdif) {
	ivertex = j;
	}
	}
	}*/

	pvertex1=pcell->vertexlist[ivertex];
	icell1=indexcell(pvertex1, pcell);
	pcellbig=pvertex1->pncell[(icell1+2)%3];


	if ((ivertexbig=indexvertex(pvertex1, pcellbig)) < 0){
		printf("In kill_cell ivertexbig < 0\n");
		exit(1);
	};

	pvertexa=pvertex1->pneighb[(icell1+1)%3];
	icella=indexcell(pvertexa, pcell);
	pvertexb=pvertex1->pneighb[(icell1+2)%3];
	icellb=indexcell(pvertexb, pcellbig);
	pvertex2=pvertex1->pneighb[icell1];
	icell2=indexcell(pvertex2, pcell);

	//check if either vertex 1 or 2 consists of only three edges
	/*
	   if (pcell->nb_vertices == 3) {
	   printf("- cell has 3 edges -> not killed\n");
	   for (i=0;i<3;i++){
	   ptmpvert=pcell->vertexlist[i];
	   for (j=0;j<3;j++){
	   pcelltmp = ptmpvert->pncell[j];
	   if (pcelltmp == pcell){
	   }
	   else if (pcelltmp->nb_vertices == 3){
	   printf("- cell has %i edges -> not killed\n",pcelltmp->nb_vertices);
	   return 0;
	   }
	   }

	   }
	   }
	   */
	for (i=0;i<3;i++){
		pcelltmp=pvertex1->pncell[(icell1+i)%3];
		if (pcelltmp->nb_vertices <=3 && pcelltmp != pcell && pcelltmp != pcellbig) {
			printf("- cell has %i edges -> not killed\n",pcelltmp->nb_vertices);
			return 0;
		}
	}
	for (i=0;i<3;i++){
		pcelltmp=pvertex2->pncell[(icell1+i)%3];
		if (pcelltmp->nb_vertices <=3 && pcelltmp != pcell && pcelltmp != pcellbig) {
			printf("- cell has %i edges -> not killed\n",pcelltmp->nb_vertices);
			return 0;
		}
	}

	pvertexc=pvertex2->pneighb[(icell2+2)%3];
	icellc=indexcell(pvertexc, pcellbig);
	pvertexd=pvertex2->pneighb[icell2];
	icelld=indexcell(pvertexd, pcell);

	pvertexa->pneighb[icella]=pvertexb;
	pvertexb->pneighb[(icellb+1)%3]=pvertexa;

	pvertexc->pneighb[icellc]=pvertexd;
	pvertexd->pneighb[(icelld+1)%3]=pvertexc;

	nbvertices1=pcell->nb_vertices;
	nbvertices2=pcellbig->nb_vertices;

	for(i=0; i < nbvertices1-2; i++){
		pvert=pcell->vertexlist[(ivertex+2+i)%nbvertices1];
		icell=indexcell(pvert, pcell);
		pvert->pncell[icell]=pcellbig;
	}

	vertexlisttemp=createnodevector(nbvertices1+nbvertices2-4);

	for(i=0; i < nbvertices1-2; i++){
		vertexlisttemp[i]=pcell->vertexlist[(ivertex+2+i)%nbvertices1];
	}
	for(i=0; i < nbvertices2-2; i++){
		vertexlisttemp[nbvertices1-2+i]=pcellbig->vertexlist[(ivertexbig+1+i)%nbvertices2];
	}


	web_dual[indexcellinwebdual(pcell)]=NULL;
	freecell(pcell);

	web[indexvertexinweb(pvertex1)]=NULL;
	web[indexvertexinweb(pvertex2)]=NULL;

	free(pcellbig->vertexlist);

	pcellbig->vertexlist=vertexlisttemp;
	pcellbig->nb_vertices=nbvertices1+nbvertices2-4;

	pcellbig->area_soll = 1.; //(1+cell_area(pcellbig))*0.5;
	//printf("kill area %f\n",pcellbig->area_soll);

	delvertex(pvertex1->pncell[(icell1+1)%3], pvertex1);
	delvertex(pvertex2->pncell[(icell2+2)%3], pvertex2);
	free(pvertex1);
	free(pvertex2);

	shrink_data();

	//snprintf(filename, 256, "./lattice.eps");
	//epsfile(filename,0);

	if (clone == 1) {
		nb_clone_tot--;
	}

	return 1;
}



/* Division functions */

void division_cell(Dual *pcell, int ivertex)
{
	int icell;
	Node *pvertex;

	pvertex=pcell->vertexlist[ivertex];
	icell=indexcell(pvertex, pcell);
	division_eq(pvertex, icell);
}

void division_eq(Node *pvertex, int icell)
{
	int n, icelltemp, ivert0cell0, nbvertices, last_vert, icell1, icell2, icell3;
	int ivertex, j,count,k;
	Node *pvert, *pvertprev, *pvertex1, *pvertex2, *pvertex3, *pnewvertex1, *pnewvertex2, **pvertexlisttemp;
	Dual *pcell, *pnewcell, *pncell0m1, *pncell2m1;
	double a, a0, x, y, z;
	vector3 testpt;

	nb_vertex_tot += 2;
	nb_cell_tot++;

	web = reallocnodevector(web, nb_vertex_tot);
	web[nb_vertex_tot - 2]=allocnode();
	web[nb_vertex_tot - 1]=allocnode();
	pnewvertex1=web[nb_vertex_tot - 2];
	pnewvertex2=web[nb_vertex_tot - 1];

	web_dual= reallocdualvector(web_dual, nb_cell_tot);
	web_dual[nb_cell_tot - 1]=allocdual();
	pnewcell=web_dual[nb_cell_tot - 1];

	pcell=pvertex->pncell[icell];


	if ((ivert0cell0=indexvertex(pvertex, pcell)) < 0){
		printf("In division ivert0cell0 < 0\n");
		exit(1);
	};

	// choose cleavage plane such that areas of daughter cells are close to 0.5
	a0 = area(pcell);
	nbvertices=pcell->nb_vertices;
	a = 0.;
	for (j=1;j <= nbvertices && a<a0/2; j++) {

		pvert    =pcell -> vertexlist[(ivert0cell0+1)%nbvertices];
		pvertprev=pcell -> vertexlist[ivert0cell0%nbvertices];
		testpt=avgpos(pvertprev,pvert);
		x = testpt.x;
		y = testpt.y; 
		z = testpt.z;
		a = part_area(pcell,j,ivert0cell0,x,y,z);//kevin_tools.c partial area calc
	}
	last_vert = j-2;
	pvertex =pcell->vertexlist[(ivert0cell0 )%nbvertices];
	pvertex1=pcell->vertexlist[(ivert0cell0+1)%nbvertices];
	pvertex2=pcell->vertexlist[(ivert0cell0+last_vert)%nbvertices];
	pvertex3=pcell->vertexlist[(ivert0cell0+last_vert+1)%nbvertices];

	icell1=indexcell(pvertex1, pcell);
	icell2=indexcell(pvertex2, pcell);
	icell3=indexcell(pvertex3, pcell);

	pncell0m1=pvertex->pncell[(icell+2)%3];

	pncell2m1=pvertex2->pncell[(icell2+2)%3];

	testpt=avgpos(pvertex,pvertex1);
	pnewvertex1->x=testpt.x;
	pnewvertex1->y=testpt.y;
	pnewvertex1->z=testpt.z;
	pnewvertex1->pneighb[0]=pvertex;
	pnewvertex1->pneighb[1]=pvertex1;
	pnewvertex1->pneighb[2]=pnewvertex2;
	pnewvertex1->pncell[0]=pncell0m1;
	pnewvertex1->pncell[1]=pnewcell;
	pnewvertex1->pncell[2]=pcell;

	testpt=avgpos(pvertex2,pvertex3);
	pnewvertex2->x=testpt.x;pnewvertex2->y=testpt.y;pnewvertex2->z=testpt.z;
	//pnewvertex2->x=(pvertex2->x + pvertex3->x)/2;
	//pnewvertex2->y=(pvertex2->y + pvertex3->y)/2;
	//pnewvertex2->z=(pvertex2->z + pvertex3->z)/2;
	pnewvertex2->pneighb[0]=pvertex2;
	pnewvertex2->pneighb[1]=pvertex3;
	pnewvertex2->pneighb[2]=pnewvertex1;
	pnewvertex2->pncell[0]=pncell2m1;
	pnewvertex2->pncell[1]=pcell;
	pnewvertex2->pncell[2]=pnewcell;


	pvertex->pneighb[icell]=pnewvertex1;

	pvertex1->pneighb[(icell1+1)%3]=pnewvertex1;

	pvertex2->pneighb[icell2]=pnewvertex2;

	pvertex3->pneighb[(icell3+1)%3]=pnewvertex2;

	for (n=1 ; n < last_vert +1 ; n++){
		pvert=pcell->vertexlist[(ivert0cell0+n)%nbvertices];
		icelltemp=indexcell(pvert, pcell);
		pvert->pncell[icelltemp]=pnewcell;
	}

	pnewcell->vertexlist=createnodevector(last_vert+2);
	pnewcell->nb_vertices=last_vert+2;
	pnewcell->marker=pcell->marker;


	pnewcell->vertexlist[0]=pnewvertex1;
	for (n=1 ; n < last_vert +1 ; n++){
		pnewcell->vertexlist[n]=pcell->vertexlist[(ivert0cell0+n)%nbvertices];
	}
	pnewcell->vertexlist[last_vert +1]=pnewvertex2;


	pvertexlisttemp=createnodevector(nbvertices-last_vert+2);
	pcell->nb_vertices=nbvertices-last_vert+2;

	pvertexlisttemp[0]=pnewvertex2;
	for (n=1 ; n < nbvertices-last_vert +1 ; n++){
		pvertexlisttemp[n]=pcell->vertexlist[(ivert0cell0+last_vert+n)%nbvertices];
	}
	pvertexlisttemp[nbvertices-last_vert+1]=pnewvertex1;

	free(pcell->vertexlist);
	pcell->vertexlist=pvertexlisttemp;

	insvertex(pncell0m1, pnewvertex1 ,pvertex);

	insvertex(pncell2m1, pnewvertex2 ,pvertex2);

	pnewcell->marker.div_rate = pcell->marker.div_rate;
	pnewcell->marker.tension_index = pcell->marker.tension_index;
	pnewcell->marker.border = pcell->marker.border;
	pnewcell->area_soll0 = pcell->area_soll0;
	pnewcell->area_soll = pnewcell->area_soll0;
	pcell->area_soll = pcell->area_soll0;

	//    pnewcell->area_soll=area(pnewcell);
	//pcell->area_soll    = 1.; //area(pcell);
	//pnewcell->area_soll = 1.;
	//printf("areas %f %f\n",area(pcell),area(pcell)/(area(pcell)+area(pnewcell)));
	count = 0;
	for(k=0;k<nb_vertex_tot;k++){
		if(web[k]->border==1){count++;}
	}
	printf("areas %f %f %d\n",area(pcell),area(pcell)/(area(pcell)+area(pnewcell)),count);
}

/* Effective division : it takes into account the probability of division and call division() */
	void
effdivision()
{
	int celltodivide; 	/* index of the cell to divide in web_dual */
	int nbvertices; 	/* # of vertices in the cell, used to determine the stating point of the division */
	int ivertex; 	/* index of the starting point of the division */
	static int itime_pressure;	/* Time counting for pressure killing process */

	Dual *pcell;


	if (gsl_rng_uniform(rng) < DIVISION_RATE){
		celltodivide = (int)(gsl_rng_uniform(rng) * (nb_cell_tot-1))+1;

		pcell=web_dual[celltodivide];
		nbvertices=pcell->nb_vertices;


		/* Checking point activation if CHECKED_DIVISION != 0 */
#if CHECKED_DIVISION
		if (!(pcell->marker.checked_division)){
#if DIVISION_RATE_CHANGEABLE
			if (pcell -> marker.div_rate==1 && gsl_rng_uniform(rng) < division_rates[1]){
				pcell->marker.checked_division=1;
			}
			else if (gsl_rng_uniform(rng) < division_rates[0]){
				pcell->marker.checked_division=1;
			}
			return;
#else

			pcell->marker.checked_division=1;
			return;
#endif /* DIVISION_RATE_CHANGEABLE */
		}
#endif /* CHECKED_DIVISION */



		/* Determination of the vertex near axis division :
		 * mechanically driven if AXIS_DIVISION != 0
		 * stochastically driven otherwise */
#if AXIS_DIVISION
		ivertex=(int)(gsl_rng_uniform(rng)*nbvertices);
#else
		ivertex = (int)(gsl_rng_uniform(rng)*nbvertices);
#endif /* AXIS_DIVISION */



#if DIVISION_RATE_CHANGEABLE
		if (pcell -> marker.div_rate==1  && gsl_rng_uniform(rng) < division_rates[1]){
			division_cell(pcell, ivertex);
			itime_pressure=1;
		}
		else if (gsl_rng_uniform(rng) < division_rates[0]){
			division_cell(pcell, ivertex);
			itime_pressure=1;
		}
#else
		division_cell(pcell, ivertex);
		itime_pressure=1;
#endif
	}

	/* Time based on stochastic behavior */
	itime++;

	itime_pressure++;

	/* Kill cells which are too much stressed */
#if PRESSURE_KILL
	if (!(itime_pressure%PRESSURE_KILL_INTER)){
		pressure_kill();
	}
#endif /* PRESSURE_KILL */
}

//counts the total number of bonds and sets all the
//pointers to bonds for each vertex to NULL to make
//the next part (where i call initbonds()) simpler
void countbonds(){
	int i,j,nb_bonds=0;
	Node *pvert;

	for(i=0;i<nb_vertex_tot;i++){
		for(j=0;j<3;j++){
			pvert = web[i]->pneighb[j];
			if(pvert->idx > i){nb_bonds++;}
		}
	}
	nb_bond_tot = nb_bonds;
}

//initializes the memory for each bond as well as
//the pointers for the endpoint vertices as well
//as bond pointers for the vertices
void initbonds(){
	int i,j,k=0,l;
	Node *pvert,*pnvert;
	Bond *pbond;

	for(i=0;i<nb_vertex_tot;i++){
		pvert = web[i];
		pvert->idx = i;
		for(j=0;j<3;j++){
			pnvert = pvert->pneighb[j];
			//sees if the neighboring vertex has greater overall index
			//if so, we create a new bond allocation and initialize
			//all the parameters as well as if the neighboring vertex
			//has index j, then the corresponding bond between the vertices
			//also has index j as a neighbor
			if(pnvert->idx > i){
				web_bond[k] = allocbond();
				pbond = web_bond[k];
				pbond->idx = k;
				pbond->pnvert[0] = pvert;
				pbond->pnvert[1] = pnvert;
				pvert->pnbond[j] = pbond;
				for(l=0;l<3;l++){
					//keeps a consistent index between neighboring vertices
					//and neighboring bonds, but in reverse as a neighbor
					if(pnvert->pneighb[l]==pvert){pnvert->pnbond[l] = pbond;}
				}
				k++;//increment the web_bond number
			}
		}
	}
	if(k!=nb_bond_tot){locerror("initbonds","Error in bond counting");}
}

//initializes the celllist and bondlist
//pointer arrays in the node structure
//also contains the algorithm for finding
//the bounding bonds of a cell
void dualneighb(){
	int i,j,k,nb_vert;
	Dual *pcell,*pcell1;
	Node *pvert,*pnvert;
	Bond *pbond;

	for(i=0;i<nb_cell_tot;i++){
		pcell = web_dual[i];
		pcell->idx = i;
		nb_vert = pcell->nb_vertices;
		pcell->celllist = alloccelllist(nb_vert);
		pcell->bondlist = allocbondlist(nb_vert);
		for(j=0;j<nb_vert;j++){
			//looks at the vertices of a cell
			pvert = pcell->vertexlist[j];
			pnvert = pcell->vertexlist[(j+1)%nb_vert];
			for(k=0;k<3;k++){
				//sees which neighbor of the vertex qualifies as the corresponding
				//neighbor in the cell.  then calls the bond between vertices j and j+1
				//the bordering bond j of the cell. can be done since above the corresp
				//bond and vertex neighbors have the same index.
				if(pnvert == pvert->pneighb[k]){
					pcell->bondlist[j] = pvert->pnbond[k];
				}
			}
		}
	}
}

//algorithm for determining the neighboring bonds and cells to a bond
void bondneighb(){
	int i,j,k,l=0,m=2,n;
	Node *pvert1, *pvert2;
	Dual *pcell1, *pcell2, *pcell3, *pcell4;
	Bond *pbond;

	for(i=0;i<nb_bond_tot;i++){
		pbond = web_bond[i];
		pvert1 = pbond->pnvert[0];pvert2 = pbond->pnvert[1];
		l=0;m=2;//initialize the array counting
		for(k=0;k<3;k++){ //this does the neighboring bonds
			if(pvert1->pnbond[k]!=pbond){pbond->pnbond[l] = pvert1->pnbond[k];l++;}
			if(pvert2->pnbond[k]!=pbond){pbond->pnbond[m] = pvert2->pnbond[k];m++;}
			if(l>2){locerror("bondneighb(l)","Error in neighboring bond assignment");}
			if(m>4){printf("m = %d, ",m);locerror("bondneighb(m)","Error in neighboring bond assignment");}
		}
		l=0; //initialize the array index
		for(k=0;k<3;k++){ //this does the neighboring cells
			n=0;
			while(n<3){ //if the neighbors match, we assign them into pncell[0] or pncell[2]
				if(pvert1->pncell[k]==pvert2->pncell[n]){pbond->pncell[l]=pvert1->pncell[k];l+=2;break;}
				n++;
			}
		}
		for(k=0;k<3;k++){ //if the neighbors aren't either of the two above, they are the ends
			if((pvert1->pncell[k]!=pbond->pncell[0]) && (pvert1->pncell[k]!=pbond->pncell[2])){
				pbond->pncell[1] = pvert1->pncell[k];
			}
			if((pvert2->pncell[k]!=pbond->pncell[0]) && (pvert2->pncell[k]!=pbond->pncell[2])){
				pbond->pncell[3] = pvert2->pncell[k];
			}
		}
	}
}

//this is here since this is much easier to do after the bond
//neighbors have been established due to our configuration
void neighbcell(){
	int i,j,k,nb_bonds;
	Bond *pbond;
	Dual *pcell,*pncell;

	for(i=0;i<nb_cell_tot;i++){
		pcell = web_dual[i];
		nb_bonds = pcell->nb_vertices;
		for(j=0;j<nb_bonds;j++){
			pbond = pcell->bondlist[j];
			if(pbond->pncell[0]!=pcell){pcell->celllist[j] = pbond->pncell[0];}
			if(pbond->pncell[2]!=pcell){pcell->celllist[j] = pbond->pncell[2];}
		}
	}
}

void createcellbondneighbs(){
	countbonds();
	web_bond = createbondvector(nb_bond_tot);
	initbonds(); //initialize bonds
	dualneighb(); //find neighbor info for cells
	bondneighb(); //find neighboring info for bonds
	neighbcell(); //find neighboring cells to cells
	if(TOROIDAL==1){toruspositioning();} //find where vertices should be located on torus
}

void updatedata(){updatecellbondneighbs();} //legacy

void updatecellbondneighbs(){
	if(web_bond!=NULL) cleanup();
	countbonds();
	web_bond = createbondvector(nb_bond_tot);
	initbonds(); //initialize bonds
	dualneighb(); //find neighbor info for cells
	bondneighb(); //find neighboring info for bonds
	neighbcell(); //find neighboring cells to cells
	if(TOROIDAL==1){toruspositioning();} //find where vertices should be located on torus
}

//cleans up certain data structures created via updatecellbondneighbs()
void cleanup(){
	int i;
	Dual *pcell;

	for(i=0;i<nb_cell_tot;i++){
		pcell = web_dual[i];
		free(pcell->bondlist);pcell->bondlist=NULL;
		free(pcell->celllist);pcell->celllist=NULL;
	}
	for(i=0;i<nb_bond_tot;i++){
		free(web_bond[i]);web_bond[i]=NULL;
	}
	free(web_bond);
	web_bond=NULL;
}

void toruspositioning(){
	int i;
	Node *pvert;

	for(i=0;i<nb_vertex_tot;i++){
		pvert=web[i];
		while(pvert->x > LX){pvert->x -= LX;}
		while(pvert->x < LXMIN){pvert->x +=LX;}
		while(pvert->y > LY){pvert->y -= LY;}
		while(pvert->y < LYMIN){pvert->y +=LY;}
	}
}

/* pressure_kill() :
 * kill cells whose pressure is too high
 */
	void
pressure_kill()
{
	int icell;	/* Loop variable */

	Dual *pcell;

	double pressure_cell;

	for (icell=1 ; icell<nb_cell_tot ; icell++){
		pcell=web_dual[icell];

		pressure_cell=pressure(pcell);

		if (pressure_cell > PRESSURE_THRESHOLD){
			kill_cell(pcell, 0);
			nb_cell_killed++;
			printf("KILLED cell %i pressure %f total killed %i\n",icell, pressure_cell, nb_cell_killed);
		}
	}
}



/* Implementation of the t1 process
 * It has to be noticed that the choice of the vertex c is not arbitrary 
 * and can induce a rotation on the lattice. This effect should be negligible though (it is in practice)
 */
	void
t1transform(Node *pvertex, int neighbor)
{
	int indexbneighbc, indexcneighbb, indexdneighba;
	Node *pvertexb, *pvertexc, *pvertexd, *pvertexe, *pvertexf; /* The explanation of the variables is on the figure in the documentation */
	Dual *pcellalpha, *pcellbeta, *pcellgamma, *pcelldelta;

	pvertexb=pvertex->pneighb[neighbor];

	pcellalpha=pvertex->pncell[neighbor];
	pcellbeta=pvertex->pncell[(neighbor+2)%3];
	indexbneighbc=indexcell(pvertexb, pcellalpha);
	pvertexc=pvertexb->pneighb[indexbneighbc];
	pvertexd=pvertex->pneighb[(neighbor+2)%3];

	pcellgamma=pvertex->pncell[(neighbor+1)%3];
	pcelldelta=pvertexb->pncell[(indexbneighbc+2)%3];

	indexcneighbb=indexcell(pvertexc, pcelldelta);
	indexdneighba=indexcell(pvertexd, pcellgamma);

	pvertexe=pvertex->pneighb[(neighbor+1)%3];
	pvertexf=pvertexb->pneighb[(indexbneighbc+2)%3];

	if(pcellalpha->nb_vertices<4 || pcellbeta->nb_vertices<4){return;}

	/* Changes */

	pvertexd->pneighb[indexdneighba]=pvertexb;
	pvertexc->pneighb[indexcneighbb]=pvertex;

	insvertex(pcellgamma, pvertexb, pvertex);
	insvertex(pcelldelta, pvertex, pvertexb);

	delvertex(pcellalpha, pvertexb);
	delvertex(pcellbeta, pvertex);

	pvertex->pneighb[0]=pvertexc;
	pvertex->pneighb[1]=pvertexe;
	pvertex->pneighb[2]=pvertexb;
	pvertex->pncell[0]=pcellalpha;
	pvertex->pncell[1]=pcellgamma;
	pvertex->pncell[2]=pcelldelta;

	pvertexb->pneighb[0]=pvertexd;
	pvertexb->pneighb[1]=pvertexf;
	pvertexb->pneighb[2]=pvertex;
	pvertexb->pncell[0]=pcellbeta;
	pvertexb->pncell[1]=pcelldelta;
	pvertexb->pncell[2]=pcellgamma;

	it1process++;
	updatecellbondneighbs();
}

void t1flagcheck(int vertexindex, int neighbindex){
	int flagt1vertex_temp, flagt1neighb_temp;
	flagt1vertex_temp = vertexindex;
	flagt1neighb_temp = neighbindex;
	if(flagt1vertex_temp == flagt1vertex && flagt1neighb_temp == flagt1neighb){flagt1++;}
	else{flagt1 = 0;}
	flagt1vertex = flagt1vertex_temp; flagt1neighb = flagt1neighb_temp;
}
