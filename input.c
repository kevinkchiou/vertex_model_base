#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "lattice.h"
#include "measurements.h"
#include "save.h"
#include "locerror.h"
#include "in_development.h"
#include "input.h"

void generatedata(FILE *file1, FILE *file2){
   int i,j,n,num;
   int buffer;
   int nbcell1,nbcell2,nbcell3;
   int nbbond1,nbbond2,nbbond3;
   int nbvert,nbvert1,nbvert2,nbvert3;
   double x,y,z,fx,fy,fz,fbuffer;
   Node *pvert;
   Dual *pcell;

   for(i=0;i<nb_vertex_tot;i++){web[i]=allocnode();}
   for(i=0;i<nb_cell_tot;i++){web_dual[i]=allocdual();}


   for(i=0;i<nb_vertex_tot;i++){
      fscanf(file1,"%d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf",&n,&nbcell1,&nbcell2,&nbcell3,&nbbond1,&nbbond2,&nbbond3,&nbvert1,&nbvert2,&nbvert3,&x,&y,&z,&fx,&fy,&fz);
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
      fscanf(file2,"%d",&n);
      fscanf(file2,"%d",&num);
      //web_dual[n] = allocdual();
      pcell = web_dual[n];
      pcell->nb_vertices = num;
      pcell->vertexlist = allocvertexlist(num);
      for(j=0;j<num;j++){
         fscanf(file2,"%d",&buffer); //don't keep the useless data
         fscanf(file2,"%d",&buffer); //don't keep the useless data
      }
      for(j=0;j<num;j++){
         fscanf(file2,"%d",&nbvert);
         pcell->vertexlist[j] = web[nbvert];
      }
      fscanf(file2,"%lf",&fbuffer);//don't keep the pressure term
      fscanf(file2,"%lf",&fbuffer);//don't keep the area term for now
      fscanf(file2,"%lf",&fbuffer);//don't keep the Asoll term for now
      fscanf(file2,"%lf",&fbuffer);//we calculate centroid positions below
      fscanf(file2,"%lf",&fbuffer);//we calculate centroid positions below
      fscanf(file2,"%lf",&fbuffer);//we calculate centroid positions below
      //populate the dual structure with data
      if(n!=0){pcell->area = area(pcell);pcell->area_soll=1.0;}
      if(n==0){pcell->area = 0.0;pcell->area_soll=0.0;}
      pcell->sqrtarea = sqrt(pcell->area);
      pcell->perimeter = perimeter(pcell);
      pcell->centroid = centroid(pcell);
      //probably might need more fscanfs in the future for the marker structure

   }
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

void infoinput(char *vertfilename, char *cellfilename){
   FILE *vertfile,*cellfile;

   //clear();
   vertfile=fopen(vertfilename,"r");cellfile=fopen(cellfilename,"r");
   if(!vertfile || !cellfile){locerror("infoinput","Cannot open files for reading!");}

   //find out the size of the arrays
   fscanf(vertfile,"%d",&nb_vertex_tot);
   fscanf(cellfile,"%d",&nb_cell_tot);

   printf("Total number of vertices = %d\n",nb_vertex_tot);
   printf("Total number of cells = %d\n",nb_cell_tot);
   //allocate the space
   web = createnodevector(nb_vertex_tot);
   web_dual = createdualvector(nb_cell_tot);
   //generate the data
   generatedata(vertfile,cellfile);
   //close the file pointers
   fclose(vertfile);fclose(cellfile);


   shrink_data();
   //test function
   //datadump();
}
