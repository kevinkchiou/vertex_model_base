#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "save.h"
#include "read.h"
#include "definelat.h"
#include "locerror.h"
#include "lattice.h"
#include "measurements.h"

int
readweb()
{
    int n, m;
    int ivert;
    char filename[256];
    FILE *pfile;
	int ct; //count for scanf and fread functions

    Node *pvertex;
    Dual *pcell;

    Nodeint nodetemp;
    Dualint dualtemp;

    printf("Enter file name : ");
    if(scanf("%s", filename)!=1){locerror("readweb()","Did not understand your entry!");}

    errno=0;

    if ( !(pfile=fopen(filename, "rb")) ){
	fprintf(stderr,"function save : Cannot open file\n lattice not saved : ");
	perror(NULL);
	return 1;
    }


    /* Number of cells and vertices */

    ct=fread(&nb_vertex_tot, sizeof(int), 1, pfile);
    ct=fread(&nb_cell_tot, sizeof(int), 1, pfile);

    web=createnodevector(nb_vertex_tot);
    web_dual=createdualvector(nb_cell_tot);

    for (n=0 ; n < nb_vertex_tot ; n++){
	web[n]=allocnode();
    }

    for (n=0 ; n < nb_cell_tot ; n++){
	web_dual[n]=allocdual();
    }

    for (n=0 ; n < nb_vertex_tot ; n++){

	pvertex=web[n];

	ct=fread(&nodetemp, sizeof(Nodeint), 1, pfile);
	pvertex->x=nodetemp.x;
	pvertex->y=nodetemp.y;
	pvertex->z=0.0;

	pvertex->pneighb[0]=web[nodetemp.ineighb[0]];
	pvertex->pneighb[1]=web[nodetemp.ineighb[1]];
	pvertex->pneighb[2]=web[nodetemp.ineighb[2]];

	pvertex->pncell[0]=web_dual[nodetemp.incell[0]];
	pvertex->pncell[1]=web_dual[nodetemp.incell[1]];
	pvertex->pncell[2]=web_dual[nodetemp.incell[2]];
    }

    for (n=0 ; n < nb_cell_tot ; n++){

	ct=fread(&dualtemp, sizeof(Dualint), 1, pfile);

	pcell=web_dual[n];

	pcell->nb_vertices=dualtemp.nb_vertices;
	pcell->marker=dualtemp.marker;

	pcell->vertexlist=createnodevector(dualtemp.nb_vertices);

	for (m=0 ; m < dualtemp.nb_vertices ; m++){
	    ct=fread(&ivert, sizeof(int), 1, pfile);
	    pcell->vertexlist[m]=web[ivert];
	}
    }

    fclose(pfile);

    return 0;
}

void generatedata(FILE *file1, FILE *file2){

	int oldstyle=0;
   int i,j,n,num;
   int buffer,clone;
   int nbcell1,nbcell2,nbcell3;
   int nbbond1,nbbond2,nbbond3;
   int nbvert,nbvert1,nbvert2,nbvert3;
   double x,y,z,fx,fy,fz,fbuffer,asoll;
   Node *pvert;
   Dual *pcell;

   int ct; //for counting number of elts read in by various scanf functions

   for(i=0;i<nb_vertex_tot;i++){web[i]=allocnode();}
   for(i=0;i<nb_cell_tot;i++){web_dual[i]=allocdual();}


   for(i=0;i<nb_vertex_tot;i++){
      ct=fscanf(file1,"%d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf",&n,&nbcell1,&nbcell2,&nbcell3,&nbbond1,&nbbond2,&nbbond3,&nbvert1,&nbvert2,&nbvert3,&x,&y,&z,&fx,&fy,&fz);
	  if(ct!=16) locerror("generatedata()","Incorrect number of elements read in!");
      //web[n] = allocnode();
      pvert = web[n];
      //populate the web structure with data
      pvert->x = x;pvert->y = y;pvert->z = z;
      pvert->pncell[0]=web_dual[nbcell1];
      pvert->pncell[1]=web_dual[nbcell2];
      pvert->pncell[2]=web_dual[nbcell3];
      pvert->pneighb[0]=web[nbvert1];
      pvert->pneighb[1]=web[nbvert2];
      pvert->pneighb[2]=web[nbvert3];
   }
   for(i=0;i<nb_vertex_tot;i++){
      pvert = web[n];
      pvert->dist_vert[0]=dist(pvert,pvert->pneighb[0]);
      pvert->dist_vert[1]=dist(pvert,pvert->pneighb[1]);
      pvert->dist_vert[2]=dist(pvert,pvert->pneighb[2]);
   }
   for(i=0;i<nb_cell_tot;i++){
      if(fscanf(file2,"%d",&n)!=1){locerror("generatedata()","Incorrect number of elements!");}
      if(fscanf(file2,"%d",&num)!=1){locerror("generatedata()","Incorrect number of elements!");}
      //web_dual[n] = allocdual();
      pcell = web_dual[n];
      pcell->nb_vertices = num;
      pcell->vertexlist = allocvertexlist(num);
      for(j=0;j<num;j++){
         ct=fscanf(file2,"%d",&buffer); //don't keep the useless data
         ct=fscanf(file2,"%d",&buffer); //don't keep the useless data
      }
      for(j=0;j<num;j++){
         ct=fscanf(file2,"%d",&nbvert);
         pcell->vertexlist[j] = web[nbvert];
      }
      ct=fscanf(file2,"%lf",&fbuffer);//don't keep the pressure term
		if(oldstyle!=1){
			ct=fscanf(file2,"%lf",&fbuffer);//don't keep the area term for now
			ct=fscanf(file2,"%lf",&asoll);//KEEP this asoll for toroid sitaution!
			ct=fscanf(file2,"%d",&clone);//KEEP clonal index
			ct=fscanf(file2,"%lf",&fbuffer);//we calculate centroid positions below
			ct=fscanf(file2,"%lf",&fbuffer);//we calculate centroid positions below
			ct=fscanf(file2,"%lf",&fbuffer);//we calculate centroid positions below
		}
		//populate the dual structure with data
		pcell->area = area(pcell);pcell->area_soll=asoll;
      pcell->marker.clone_index=clone;
      if(n==0 && TOROIDAL==0){pcell->area = 0.0;pcell->area_soll=0.0;}
      pcell->sqrtarea = sqrt(pcell->area);
      pcell->perimeter = perimeter(pcell);
      pcell->centroid = centroid(pcell);
      //probably might need more fscanfs in the future for the marker structure

   }
   updatecellbondneighbs();//to create other parts of the structure!
}

//clears the web and web_dual caused by coord_gen so i can get going
void clear(){
   int i,j,num;
   Dual *pcell;
   Node *pvert;

   for(i=0;i<nb_vertex_tot;i++){
      free(web[i]);
   }
   free(web);

   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      num = pcell->nb_vertices;
      for(j=0;i<num;i++){
         free(pcell->vertexlist[j]);
      }
      free(pcell);
   }
   free(web_dual);
}

void datadump(){
   int i,j,num;
   Dual *pcell;
   Node *pvert;
   FILE *pfile1,*pfile2;

   pfile1 = fopen("vertexdump.dat","w");pfile2=fopen("celldump.dat","w");

   for(i=0;i<nb_vertex_tot;i++){
      pvert = web[i];
      fprintf(pfile1,"%d\t",i);
      for(j=0;j<3;j++){
         if(j!=2){fprintf(pfile1,"%d ",indexcellinwebdual(pvert->pncell[j]));}
         if(j==2){fprintf(pfile1,"%d\t",indexcellinwebdual(pvert->pncell[j]));}
      }
      for(j=0;j<3;j++){
         if(j!=2){fprintf(pfile1,"%d ",indexvertexinweb(pvert->pneighb[j]));}
         if(j==2){fprintf(pfile1,"%d\t",indexvertexinweb(pvert->pneighb[j]));}
      }
      fprintf(pfile1,"\n");
   }
   fclose(pfile1);

   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      num = pcell->nb_vertices;
      fprintf(pfile2,"%d\t%d\t",i,num);
      for(j=0;j<num;j++){
         if(j!=num-1){fprintf(pfile2,"%d ",indexvertexinweb(pcell->vertexlist[j]));}
         if(j==num-1){fprintf(pfile2,"%d\t",indexvertexinweb(pcell->vertexlist[j]));}
      }
      fprintf(pfile2,"\n");
   }
   fclose(pfile2);

}

void infoinput(char *vertfilename, char *cellfilename,int torus){
   FILE *vertfile,*cellfile;

   TOROIDAL=torus;
   if(TOROIDAL==1){
       LX=20.;LY=LX*sqrt(3)*0.5;
       LXMIN=0.;LYMIN=0.;
   }
   if(TOROIDAL==2){
      TOROIDAL=1;
      LX=40.;LY=LX*sqrt(3)*0.25;
      LXMIN=0.;LYMIN=0.;
   }
   //clear();
   vertfile=fopen(vertfilename,"r");cellfile=fopen(cellfilename,"r");
   if(!vertfile || !cellfile){locerror("infoinput","Cannot open files for reading!");}

   //find out the size of the arrays
   if(fscanf(vertfile,"%d",&nb_vertex_tot)!=1){locerror("infoinput()","Incorrect number of elements scanned in!");}
   if(fscanf(cellfile,"%d",&nb_cell_tot)!=1){locerror("infoinput()","Incorrect number of elements scanned in!");}

   printf("Total number of vertices = %d\n",nb_vertex_tot);
   printf("Total number of cells = %d\n",nb_cell_tot);
   //allocate the space
   web = createnodevector(nb_vertex_tot);
   web_dual = createdualvector(nb_cell_tot);
   //generate the data
   generatedata(vertfile,cellfile);
   //close the file pointers
   fclose(vertfile);fclose(cellfile);


   //test function
   //datadump();
}
