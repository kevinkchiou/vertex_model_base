
#include <stdio.h>
#include <errno.h>

#include "definelat.h"
#include "lattice.h"
#include "save.h"
#include "read.h"
#include "measurements.h"
#include "locerror.h"
//#include "kevin_tools.h"


/* openfile function
 * It eases the error mechanism
 */
FILE
*openfile(char *filename)
{
    FILE *pfile;

    errno=0;

    if (!(pfile=fopen(filename, "w"))){
	fprintf(stderr,"function openfile : Cannot open file\n");
	perror(NULL);
	return NULL;
    }

    return pfile;
}


/* Record the lattice information in a binary file
 */
int
saveweb(char *filename)
{
    int n, m;
    int ivert;
    FILE *pfile;

    Node *pvertex;
    Dual *pcell;

    Nodeint nodetemp;
    Dualint dualtemp;


    pfile=openfile(filename);
    if (!pfile){
	fprintf(stderr,"Cannot save lattice\n");
	return 1;
    }


    /* Number of cells and vertices */

    fwrite(&nb_vertex_tot, sizeof(int), 1, pfile);
    fwrite(&nb_cell_tot, sizeof(int), 1, pfile);

    for (n=0 ; n < nb_vertex_tot ; n++){

	pvertex=web[n];

	nodetemp.x=pvertex->x;
	nodetemp.y=pvertex->y;

	nodetemp.ineighb[0]=indexvertexinweb(pvertex->pneighb[0]);
	nodetemp.ineighb[1]=indexvertexinweb(pvertex->pneighb[1]);
	nodetemp.ineighb[2]=indexvertexinweb(pvertex->pneighb[2]);

	nodetemp.incell[0]=indexcellinwebdual(pvertex->pncell[0]);
	nodetemp.incell[1]=indexcellinwebdual(pvertex->pncell[1]);
	nodetemp.incell[2]=indexcellinwebdual(pvertex->pncell[2]);

	fwrite(&nodetemp, sizeof(Nodeint), 1, pfile);
    }

    for (n=0 ; n < nb_cell_tot ; n++){

	pcell=web_dual[n];

	dualtemp.nb_vertices=pcell->nb_vertices;
	dualtemp.marker=pcell->marker;

	fwrite(&dualtemp, sizeof(Dualint), 1, pfile);

	for (m=0 ; m < dualtemp.nb_vertices ; m++){
	    ivert=indexvertexinweb(pcell->vertexlist[m]);
	    fwrite(&ivert, sizeof(int), 1, pfile);
	}
    }

    fclose(pfile);

    return 0;
}

void latticeprint(char *filename, unsigned int type){
	int i,ivert;
	FILE *pfile;
	Dual *pcell;
	Node *pvert;
	vector3 center,centroid;
	double xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
	
	pfile = openfile(filename);
	fprintf(pfile,"\n");
	
	for(i=1;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		fprintf(pfile,"%lf\n",pressure(pcell));
		pvert=pcell->vertexlist[0];
		for(ivert=0;ivert<pcell->nb_vertices;ivert++){
			pvert=pcell->vertexlist[ivert];
			if(ivert==0){
				xmin = pvert->x;
				xmax = pvert->x;
				ymin = pvert->y;
				ymax = pvert->y;
				zmin = pvert->z;
				zmax = pvert->z;
			}
			else{
				if (xmin > pvert->x) xmin = pvert->x;
				if (xmax < pvert->x) xmax = pvert->x;
				if (ymin > pvert->y) ymin = pvert->y;
				if (ymax < pvert->y) ymax = pvert->y;
				if (zmin > pvert->z) zmin = pvert->z;
				if (zmax < pvert->z) zmax = pvert->z;
			}
			x=pvert->x;y=pvert->y;z=pvert->z;
			fprintf(pfile,"%lf %lf %lf\n",x,y,z);
		}
		fprintf(pfile,"\n");
		fprintf(pfile,"\n");
	}
	fprintf(pfile,"END\n");
	fclose(pfile);
}

void infoprint(char *filename, unsigned int type){
   int i,j,k,num;
   char filename1[256],filename2[256],filename3[256];
   FILE *pfile1,*pfile2,*pfile3;
   Dual *pcell,*pcell0,*pcell2;
   Node *pvert;
   Bond *pbond;
   vector3 center;
   double xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
   double test,stress_tr,stress_det;
   int *list;
   matrix22 stress;

   snprintf(filename1,255,"%s.cellinfo",filename);
   snprintf(filename2,255,"%s.bondinfo",filename);
   snprintf(filename3,255,"%s.vertinfo",filename);

   pfile1=openfile(filename1);pfile2=openfile(filename2);pfile3=openfile(filename3);
   fprintf(pfile1,"%d\n",nb_cell_tot);
   fprintf(pfile2,"%d\n",nb_bond_tot);
   fprintf(pfile3,"%d\n",nb_vertex_tot);
   for(i=0;i<nb_cell_tot;i++){
      pcell=web_dual[i];
      num = pcell->nb_vertices;
      fprintf(pfile1,"%d\t%d\t",i,num);
      for(j=0;j<num;j++){
         if(j!=num-1){fprintf(pfile1,"%d ",indexcellinwebdual(pcell->celllist[j]));}
         if(j==num-1){fprintf(pfile1,"%d\t",indexcellinwebdual(pcell->celllist[j]));}
      }
      for(j=0;j<num;j++){
         if(j!=num-1){fprintf(pfile1,"%d ",indexbondinwebbond(pcell->bondlist[j]));}
         if(j==num-1){fprintf(pfile1,"%d\t",indexbondinwebbond(pcell->bondlist[j]));}
      }
      for(j=0;j<num;j++){
         if(j!=num-1){fprintf(pfile1,"%d ",indexvertexinweb(pcell->vertexlist[j]));}
         if(j==num-1){fprintf(pfile1,"%d\t",indexvertexinweb(pcell->vertexlist[j]));}
      }
      center=centroid(pcell);
      stress=stress_cells_2d(&i,1);
	  stress_tr=stress.c[0][0]+stress.c[1][1];
	  stress_det=stress.c[0][0]*stress.c[1][1]-stress.c[0][1]*stress.c[1][0];
      fprintf(pfile1,"%lf %lf %lf %lf %lf %d ",stress_det,stress_tr,pressure(pcell),area(pcell),pcell->area_soll,pcell->marker.border);
      fprintf(pfile1,"%lf %lf %lf\n",center.x,center.y,center.z);
   }
   fclose(pfile1);
   for(i=0;i<nb_bond_tot;i++){
      pbond = web_bond[i];
      fprintf(pfile2,"%d\t",i);
      for(j=0;j<4;j++){
         if(j!=3){fprintf(pfile2,"%d ",indexcellinwebdual(pbond->pncell[j]));}
         if(j==3){fprintf(pfile2,"%d\t",indexcellinwebdual(pbond->pncell[j]));}
      }
      for(j=0;j<4;j++){
         if(j!=3){fprintf(pfile2,"%d ",indexbondinwebbond(pbond->pnbond[j]));}
         if(j==3){fprintf(pfile2,"%d\t",indexbondinwebbond(pbond->pnbond[j]));}
      }
      for(j=0;j<2;j++){
         if(j!=1){fprintf(pfile2,"%d ",indexvertexinweb(pbond->pnvert[j]));}
         if(j==1){fprintf(pfile2,"%d\t",indexvertexinweb(pbond->pnvert[j]));}
      }
      fprintf(pfile2,"%lf\t",dist(pbond->pnvert[0],pbond->pnvert[1]));
      pcell0=pbond->pncell[0];pcell2=pbond->pncell[2];
      fprintf(pfile2,"%lf\n",tension(pbond));
   }
   fclose(pfile2);
   for(i=0;i<nb_vertex_tot;i++){
      pvert = web[i];
      fprintf(pfile3,"%d\t",i);
      for(j=0;j<3;j++){
         if(j!=2){fprintf(pfile3,"%d ",indexcellinwebdual(pvert->pncell[j]));}
         if(j==2){fprintf(pfile3,"%d\t",indexcellinwebdual(pvert->pncell[j]));}
      }
      for(j=0;j<3;j++){
         if(j!=2){fprintf(pfile3,"%d ",indexbondinwebbond(pvert->pnbond[j]));}
         if(j==2){fprintf(pfile3,"%d\t",indexbondinwebbond(pvert->pnbond[j]));}
      }
      for(j=0;j<3;j++){
         if(j!=2){fprintf(pfile3,"%d ",indexvertexinweb(pvert->pneighb[j]));}
         if(j==2){fprintf(pfile3,"%d\t",indexvertexinweb(pvert->pneighb[j]));}
      }
      fprintf(pfile3,"%.9f %.9f %.9f\t",pvert->x,pvert->y,pvert->z);
      fprintf(pfile3,"%lf %lf %lf\n",forcevector[0][i],forcevector[1][i],forcevector[2][i]);
   }
   fclose(pfile3);
}

int createFilenames(size_t n,char *svgfname,char *infofname,char *svgdir,char *infodir,char *fname){
	snprintf(svgfname,n,"%s/%s.svg",svgdir,fname);
	snprintf(infofname,n,"%s/%s",infodir,fname);
	return 0;

}
