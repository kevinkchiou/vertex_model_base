/* The comments directly come from the libpng web site...
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "definelat.h"
#include "lattice.h"
#include "measurements.h"
#include "save.h"
#include "pngwrite.h"

FILE *fp;


/* Svg file generation : takes a snapshop of the lattice
 * It is very usefull to test the program and/or to know what happen when there is a bug
 */
void
svglattice(char *filename, unsigned int type)
{
   int i,start=1;
    FILE *pfile;
    Dual *pcell;
    vector3 center;

    /* File initialisation */
    pfile = openfile(filename);

    fprintf(pfile, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n\
<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n\
   \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\
<svg\n   xmlns=\"http://www.w3.org/2000/svg\"\n\
   version=\"1.1\"\n\
   x=\"0mm\"\n\
   y=\"0mm\"\n\
   width=\"250mm\"\n\
   height=\"250mm\">\n\
   <desc>Growing lattice\n\
   </desc>\n\
   <defs>\n\
   </defs>\n\
   <g transform=\"translate(350,350)\">\n\
   <g transform=\"scale(20)\"\n\
      style=\"fill:white;stroke:black;stroke-width:0.05;stroke-linecap:round;stroke-linejoin:round\">\n");


    fprintf(pfile, "<polygon id=\"cell0\" fill = \"rgb(230,230,230)\" stroke=\"#000000\" stroke-width=\"0.01\" points =\" -12,-12 -12,12 12,12 12,-12 \" />\n");

    center.x = 0.;
    center.y = 0.;
	if(TOROIDAL>0){start=0;}
    for (i=start;i<nb_cell_tot;i++){
      pcell=web_dual[i];
      pcell->centroid = centroid(pcell);
      center.x += pcell->centroid.x;
      center.y += pcell->centroid.y;
    }
    center.x /= (nb_cell_tot-start);
    center.y /= (nb_cell_tot-start);
	if(TOROIDAL>0){center.x=LX/2.;center.y=LY/2.;}
    for (i=start;i<nb_cell_tot;i++){
       pcell=web_dual[i];
       svgprcell(pcell,center,pfile,type);
    }
    fprintf(pfile, "      </g>\n   </g>\n</svg>");
    fclose(pfile);
}


void
svgprcell(Dual *pcell, vector3 center, FILE *pfile, unsigned int type)
{
  int ivert;
  Node *pvert;
  double xmin,xmax,ymin,ymax,pres;

  pvert=pcell->vertexlist[0];

  for (ivert=1; ivert< pcell->nb_vertices; ivert++){
    if (ivert == 1) {
      xmin = pvert->x-center.x;
      xmax = pvert->x-center.x;
      ymin = pvert->y-center.y;
      ymax = pvert->y-center.y;
    }
    else {
      if (xmin > pvert->x-center.x) xmin = pvert->x-center.x;
      if (xmax < pvert->x-center.x) xmax = pvert->x-center.x;
      if (ymin > pvert->y-center.y) ymin = pvert->y-center.y;
      if (ymax < pvert->y-center.y) ymax = pvert->y-center.y;

    }
  }

  fprintf(pfile,"\n<linearGradient id=\"gracell%i\" gradientUnits=\"userSpaceOnUse\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\">\n",indexcellinwebdual(pcell),xmin,ymin,xmax,ymax);
  switch(type){
     case 0:
	printdefects(pcell, pfile);
	break;
     case 1:
	printpressure(pcell, pfile);
	break;
     case 2:
	printcellstate(pcell,pfile);
	break;
     case 3:
        break;
     case 4:
	break;
     case 5:
        break;
  }

  if (pcell->nb_vertices > 6) {
    fprintf(pfile,"<stop  offset=\"0\" style=\"stop-color:#EF7070\"/>\n");
  } 
  else if (pcell->nb_vertices < 6) {
    fprintf(pfile,"<stop  offset=\"0\" style=\"stop-color:#7070EF\"/>\n");

  }
  else {
    fprintf(pfile,"<stop  offset=\"0\" style=\"stop-color:#EFEFEF\"/>\n");

  }
  fprintf(pfile,"<stop  offset=\"1\" style=\"stop-color:#707070\"/>\n");
  fprintf(pfile,"</linearGradient>\n");

  fprintf(pfile, "<polygon id=\"cell%i\"", indexcellinwebdual(pcell));
  //  fprintf(pfile, " fill=\"url(#gracell%i)\" stroke=\"#000000\" stroke-width=\"0.5\" points =\"", indexcellinwebdual(pcell));
   pres=pressure(pcell);
   
   if (pres<0.0) fprintf(pfile, " fill=\"rgb(%1$i,%1$i,255)\" stroke=\"#000000\" stroke-width=\"0.05\" points =\"", (int)(255*(1+0.01*(pres))),indexcellinwebdual(pcell));
   else fprintf(pfile, " fill=\"rgb(255,%1$i,%1$i)\" stroke=\"#000000\" stroke-width=\"0.05\" points =\"", (int)(255*(1-0.01*(pres))),indexcellinwebdual(pcell));

  pcell->centroid=centroid(pcell);
  for (ivert=0; ivert< pcell->nb_vertices; ivert++){
    pvert=pcell->vertexlist[ivert];
    fprintf(pfile, " %f,%f", ptptdistx(pvert->x,pcell->centroid.x)-center.x+pcell->centroid.x,ptptdisty(pvert->y,pcell->centroid.y)-center.y+pcell->centroid.y);
  }
  fprintf(pfile, " \" />\n");

  if (pcell->marker.clone_index==1){
     fprintf(pfile, "            <circle cx=\"%.2f\" cy=\"%.2f\" r=\"0.1\" fill=\"red\" />",\
	   pcell->centroid.x-center.x, pcell->centroid.y-center.y);
  }
  if(pcell->marker.clone_index==2){
     fprintf(pfile, "            <circle cx=\"%.2f\" cy=\"%.2f\" r=\"0.1\" fill=\"blue\" />",\
	   pcell->centroid.x-center.x, pcell->centroid.y-center.y);
  }
}

void
printpressure(Dual *pcell, FILE *pfile)
{
   double pres;
   
   pres=pressure(pcell);
   
   if (pres<12.0) fprintf(pfile, "            fill=\"rgb(%1$i,%1$i,255)\"\n", (int)(255*(1+0.2*(pres-12.0))));
   else fprintf(pfile, "            fill=\"rgb(255,%1$i,%1$i)\"\n", (int)(255*(1-0.2*(pres-12.0))));
}

void printcellstate(Dual *pcell, FILE *pfile){

	int idx=pcell->marker.clone_index;

	if(idx==1){
		fprintf(pfile, "            fill=\"rgb(255,0,0)\"\n");
	}
	if(idx==2){
		fprintf(pfile, "            fill=\"rgb(0,0,255)\"\n");
	}
	if(idx==0){
		fprintf(pfile, "            fill=\"rgb(255,255,255)\"\n");
	}
}

void
printdefects(Dual *pcell, FILE *pfile)
{
   switch (pcell->nb_vertices){
   case 4:
   case 5:
	 fprintf(pfile, "            fill=\"rgb(255,150,150)\"\n");
	 break;
      case 6:
	 fprintf(pfile, "            fill=\"rgb(182,182,182)\"\n");
	 break;
      case 7:
      case 8:
      case 9:
      case 10:
      case 11:
      case 12:
      case 13:
      case 14:
	 fprintf(pfile, "            fill=\"rgb(200,200,255)\"\n");
	 break;
   }
}



void
print_cell_eps(Dual *pcell, FILE *pfile)
{
  int ivertex;
  int nbvertices;
  double press;
  double dist, lcol;
  vector3 centroid_cell;
  Node *pvertex;

  nbvertices=pcell->nb_vertices;

  centroid_cell=centroid(pcell);
  dist = centroid_cell.x*centroid_cell.x + centroid_cell.y*centroid_cell.y;


  pvertex=pcell->vertexlist[0];
  fprintf(pfile, "n\n");
  fprintf(pfile, "%f %f m\n", pvertex -> x, pvertex -> y);

  for (ivertex=1; ivertex < nbvertices; ivertex++){
    pvertex=pcell->vertexlist[ivertex];
    fprintf(pfile, "%f %f l\n", pvertex -> x, pvertex -> y);
  }

  fprintf(pfile, "c\n");
  fprintf(pfile, "gsave\n");

  fprintf(pfile, "stroke\n");
}

/* Eps file generation : takes a snapshop of the lattice
 * It is very usefull to test the program and/or to know what happen when there is a bug
 */
void
epsfile(char *filename, unsigned int type)
{
    FILE *pfile;
    int icell;	/* Loop variables */
    time_t epstime;

    /* File initialisation */
    pfile = openfile(filename);

    epstime=time(NULL);
    /* Header of the postscript file */
    /* Comments */
    fprintf(pfile, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(pfile, "%%%%Title: %s\n", filename);
    fprintf(pfile, "%%%%Creator: Morphogenesis\n");
    fprintf(pfile, "%%%%CreationDate: %s", asctime(gmtime(&epstime)));
    fprintf(pfile, "%%%%BoundingBox: -200 -200 300 300\n");
    fprintf(pfile, "%%%%EndComments\n");

    fprintf(pfile, "/unit {10 mul} def\n");
/*    fprintf(pfile, "100 100 translate\n");*/
    fprintf(pfile, "/clonedot {\ngsave\nnewpath\n0 setgray\nexch unit exch unit translate\n0 0 .08 unit 0 360 arc closepath fill\ngrestore\n}def\n");
    fprintf(pfile, "/Times-Roman findfont\n");
    fprintf(pfile, "/ccol{0 0 0 setrgbcolor}def\n");
    fprintf(pfile, "/lcol{0 0 0 setrgbcolor}def\n");
    fprintf(pfile, "/cslw{0.5 setlinewidth}def\n");
    fprintf(pfile, "/lslw{0.5 setlinewidth}def\n");
    fprintf(pfile, "/m{exch unit exch unit moveto}def\n");
    fprintf(pfile, "/l{exch unit exch unit lineto}def\n");
    fprintf(pfile, "/n{newpath}def\n");
    fprintf(pfile, "/c{closepath}def\n");
    fprintf(pfile, "/s{stroke}def\n");
    fprintf(pfile, "/f{fill}def\n");
    fprintf(pfile, "/pf{2.5}def\n");
    fprintf(pfile, "/scor{neg 2 mul 3 add exch neg 2 mul 3 add setrgbcolor}def\n");
    fprintf(pfile, "/scob{neg 1 add exch neg 1 add  1. setrgbcolor}def\n");

    fprintf(pfile, "6 scalefont\n");
    fprintf(pfile, "setfont\n");
    fprintf(pfile, "0.2 setlinewidth\n");
    

    /* Draw all cells */
    for (icell=1 ; icell<nb_cell_tot; icell++){
	  print_cell_eps(web_dual[icell], pfile);
    }

    /* File closure */
    fprintf(pfile, "showpage\n");
    fprintf(pfile, "%%%%Trailer\n");
    fprintf(pfile, "%%EOF\n");
    fclose(pfile);
}

