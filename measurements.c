#include<stdio.h>
#include<stdlib.h>

#include<math.h>
#include <gsl/gsl_multimin.h>

#include "const.h"
#include "definelat.h"
#include "lattice.h"
#include "measurements.h"
#include "save.h"
#include "opti.h"
#include "in_development.h"

/* distance between 2 vertices */
double dist(Node *pvertex1, Node *pvertex2)
{
   register double dx, dy, dz;

   dx=distx(pvertex1,pvertex2);
   dy=disty(pvertex1,pvertex2);
   dz=distz(pvertex1,pvertex2);
   return sqrt(dx*dx + dy*dy + dz*dz);
}

/*calculates the scalar area of a cell specified */
double area(Dual *pcell){
	int ivert,nbvert;
	vector3 av,avtemp,center;
	Node *pv,*pvprev;
	
	nbvert = pcell->nb_vertices;
	pvprev = pcell->vertexlist[nbvert-1];
	pv = pcell->vertexlist[0];
	center = centroid(pcell);
	av = area_calc(pvprev,pv,center);
	
    //printf("calculating area for cell %d...",indexcellinwebdual(pcell));fflush(stdout);
	for(ivert=1;ivert<nbvert;ivert++){
		pvprev=pv;
		pv=pcell->vertexlist[ivert];
		avtemp=area_calc(pvprev,pv,center);
		av.x+=avtemp.x;
		av.y+=avtemp.y;
		av.z+=avtemp.z;
	}
	return 0.5*sqrt(pow(av.x,2)+pow(av.y,2)+pow(av.z,2));
}
	
//calculates the partial area of a cell specified starting at ivert0cell0
//and coords (x,y,z) up to a maximum number of vertices max
double part_area(Dual *pcell, int max, int ivert0cell0, double x, double y, double z){
    int ivert,nbvert;
    vector3 av,avtemp,center,point;
    Node *pv,*pvprev;
    double x1,x2,y1,y2,z1,z2;

    nbvert = pcell->nb_vertices;
    pvprev = pcell->vertexlist[(ivert0cell0)%nbvert];
    pv = pcell->vertexlist[(ivert0cell0+1)%nbvert];
    center = centroid(pcell);

    x1=ptptdistx(x,center.x);x2=-centdistx(center,pv);
    y1=ptptdisty(y,center.y);y2=-centdisty(center,pv);
    z1=ptptdistz(z,center.z);z2=-centdistz(center,pv);
    av.x = y1*z2-z1*y2;
    av.y = z1*x2-x1*z2;
    av.z = x1*y2-y1*x2;

    for(ivert=1;ivert<max;ivert++){
        pvprev=pv;
        pv=pcell->vertexlist[(ivert0cell0+ivert+1)%nbvert];
        avtemp=area_calc(pvprev,pv,center);
        av.x+=avtemp.x;
        av.y+=avtemp.y;
        av.z+=avtemp.z;
    }
    pvprev=pv;
    pv=pcell->vertexlist[(ivert0cell0+max)%nbvert];
    point=avgpos(pvprev,pv);//finds the average position between the two nodes
    x=point.x;y=point.y;z=point.z;
    //note here we switched the order to maintain orientation
    x2=ptptdistx(x,center.x);x1=-centdistx(center,pv);
    y2=ptptdisty(y,center.y);y1=-centdisty(center,pv);
    z2=ptptdistz(z,center.z);z1=-centdistz(center,pv);
    av.x += y1*z2-z1*y2;
    av.y += z1*x2-x1*z2;
    av.z += x1*y2-y1*x2;
    return 0.5*sqrt(pow(av.x,2)+pow(av.y,2)+pow(av.z,2));
}

/*calculates area vector of triangle spanned by points represented
  by the triangle calculated by the cross product pv1 x pv2 about
  the point centroid */
vector3 area_calc(Node *pv1, Node *pv2, vector3 centroid){
    double x1,x2,y1,y2,z1,z2;
    vector3 av; //area vector

    x1=-centdistx(centroid,pv1);
    x2=-centdistx(centroid,pv2);
    y1=-centdisty(centroid,pv1);
    y2=-centdisty(centroid,pv2);
    z1=-centdistz(centroid,pv1);
    z2=-centdistz(centroid,pv2);

    av.x = y1*z2-z1*y2;
    av.y = z1*x2-x1*z2;
    av.z = x1*y2-y1*x2;

    return av;
}

double distx(Node *pv1, Node *pv2){
    register double x1,x2,dx;

    x1=pv1->x;x2=pv2->x;
    dx=ptptdistx(x1,x2);

    return dx;
}
double disty(Node *pv1, Node *pv2){
    register double dy,y1,y2;

    y1=pv1->y;y2=pv2->y;
    dy=ptptdisty(y1,y2);

    return dy;

}
double distz(Node *pv1, Node *pv2){
    register double dz,z1,z2;

    z1=pv1->z;z2=pv2->z;
    dz=ptptdistz(z1,z2);
    return dz;

}
double centdistx(vector3 cent,Node *pv){
	register double dx,x1,x2;
	
	x1=cent.x;x2=pv->x;
    dx=ptptdistx(x1,x2);
	return dx;
}
double centdisty(vector3 cent,Node *pv){
	register double dy,y1,y2;
	
	y1=cent.y;y2=pv->y;
    dy=ptptdisty(y1,y2);
	return dy;
}
double centdistz(vector3 cent,Node *pv){
	register double dz,z1,z2;
	
	z1=cent.z;z2=pv->z;
    dz=ptptdistz(z1,z2);
	return dz;
}
double ptptdistx(double x1,double x2){
    register double dx,dx1,dx2,dxsqr,dx1sqr,dx2sqr;

    dx=x1-x2;
    if(TOROIDAL!=0){
        dx1=dx+LX;dx2=dx-LX;
        dxsqr=dx*dx;dx1sqr=dx1*dx1;dx2sqr=dx2*dx2;
        if(dx1sqr<dxsqr){dx=dx1;}
        if(dx2sqr<dxsqr){dx=dx2;}
    }
    return dx;
}
double ptptdisty(double y1,double y2){
    register double dy,dy1,dy2,dysqr,dy1sqr,dy2sqr;

    dy=y1-y2;
    if(TOROIDAL!=0){
        dy1=dy+LY;dy2=dy-LY;
        dysqr=dy*dy;dy1sqr=dy1*dy1;dy2sqr=dy2*dy2;
        if(dy1sqr<dysqr){dy=dy1;}
        if(dy2sqr<dysqr){dy=dy2;}
    }
    return dy;
}
double ptptdistz(double z1,double z2){
    double dz;

    dz=z1-z2;
    return dz;
}
double fixxloc(double x,double xbound){
	while(x<0.){x+=xbound;}
	while(x>xbound){x-=xbound;}

    return x;
}
double fixyloc(double y,double ybound){
	while(y<0.){y+=ybound;}
	while(y>ybound){y-=ybound;}

    return y;
}

vector3 avgpos(Node *pv1, Node *pv2){
    vector3 r;
    double dx,dy,dxtest,dytest;

    r.x=0.5*(pv1->x + pv2->x);
    r.y=0.5*(pv1->y + pv2->y);
    r.z=0.5*(pv1->z + pv2->z);

    if(TOROIDAL!=0){
        dx=distx(pv1,pv2);dy=disty(pv1,pv2);
        dxtest=centdistx(r,pv1);dytest=centdisty(r,pv1);
        if(dxtest*dxtest>dx*dx){r.x+=LX*0.5;}
        if(dytest*dytest>dy*dy){r.y+=LY*0.5;}
        while(r.x<0.){r.x+=LX;}
        while(r.x>LX){r.x-=LX;}
        while(r.y<0.){r.y+=LY;}
        while(r.y>LY){r.y-=LY;}
    }

    return r;

}

vector3 centroid(Dual *pcell){
	int i,nbneighb,count,countx,county;
	Node *pv;
	vector3 center,centertest,result;
	double tempx,tempy;
	double sumtestx,sumtesty,sumresx,sumresy;

	center.x=0;center.y=0;center.z=0;
	nbneighb=pcell->nb_vertices;

	for(i=0;i<nbneighb;i++){
		pv=pcell->vertexlist[i];

		center.x+=pv->x;
		center.y+=pv->y;
		center.z+=pv->z;
	}

	center.x /= nbneighb;
	center.y /= nbneighb;
	center.z /= nbneighb;

	if(TOROIDAL!=0){
		countx=0;county=0; //initially there's no shift
		count=-ceil(nbneighb/2);result=center; //initial things to go by.
		while(count<=ceil(nbneighb/2)){
			centertest.x=center.x+count*LX/nbneighb;
			centertest.y=center.y+count*LY/nbneighb;
			sumtestx=0;sumtesty=0;sumresx=0;sumresy=0;
			for(i=0;i<nbneighb;i++){
				pv=pcell->vertexlist[i];
				tempx=centdistx(centertest,pv);tempy=centdisty(centertest,pv);
				sumtestx+=tempx*tempx;
				sumtesty+=tempy*tempy;
				tempx=centdistx(result,pv);tempy=centdisty(result,pv);
				sumresx+=tempx*tempx;
				sumresy+=tempy*tempy;
			}
			if(sumtestx<sumresx){result.x=centertest.x;countx=count;}//find shift with smallest x-distance sum
			if(sumtesty<sumresy){result.y=centertest.y;county=count;}//find shift with smallest y-distance sum
			count++;
		}
        while(result.x>LX){result.x-=LX;}
        while(result.x<0.){result.x+=LX;}
        while(result.y>LY){result.y-=LY;}
        while(result.y<0.){result.y+=LY;}
		center=result;
	}
	return center;
}


/* Computes the pressure of a cell
 */
double pressure(Dual *pcell)
{
   return -1.0*area_fctp(pcell);
}

double tension(Bond *pbond){
	
	double t;

	//perimeter contributions
	t = rho_fctpp(pbond->pncell[0])+rho_fctpp(pbond->pncell[2]);
	//bond contributions
	t+= rhol_fctpl(pbond);
	return t;
}

/* Computes the pressure of a cell*/

double pressureavg()
{
   int i;
   double avg=0;
   Dual *pcell;
   update_pa();
   for (i=1;i<nb_cell_tot;i++){
      pcell=web_dual[i];
      avg+=pcell->area*pressure(pcell);
   }
   avg/=nb_cell_tot;
   return avg;
}


/* Computes the perimeter of a cell */

double perimeter(Dual *pcell)
{
    int i, nbvertices;
    double perim=0;


    nbvertices= pcell -> nb_vertices;

    for (i=0; i<nbvertices-1; i++){
	perim += dist(pcell -> vertexlist[i], pcell -> vertexlist[i+1]);
    }
   
    perim += dist(pcell -> vertexlist[0], pcell -> vertexlist[nbvertices-1]);

    return perim;
}

 
vector3 center_lattice()
{
    int icell;	// Loop variable
    double left, right, bottom, top, up, down;
    vector3 center, centroid_cell;
    Dual *pcell;

    pcell=web_dual[1];

    centroid_cell=centroid(pcell);
    left=right=centroid_cell.x;
    bottom=top=centroid_cell.y;
	up=down=centroid_cell.z;

    for (icell=2 ; icell < nb_cell_tot ; icell++){
	centroid_cell=centroid(web_dual[icell]);
	if (centroid_cell.x < left) left=centroid_cell.x;
	if (centroid_cell.x > right) right=centroid_cell.x;
	if (centroid_cell.y < bottom) bottom=centroid_cell.y;
	if (centroid_cell.y > top) top=centroid_cell.y;
	if (centroid_cell.z < down) bottom=centroid_cell.z;
	if (centroid_cell.z > up) top=centroid_cell.z;
    }

    center.x=0.5*(left+right);
    center.y=0.5*(bottom+top);
	center.z=0.5*(down+up);

    return center;
}

matrix22 stress_cells_2d(int *list,int length){

	int i,j,nbneighb;
	Dual *pcell;
	Bond *pbond;
	double div,p,t,dx,dy,d;
	matrix22 sig;

	div=0.;
	sig.c[0][0]=0.;sig.c[0][1]=0.;sig.c[1][0]=0.;sig.c[1][1]=0.;

	for(i=0;i<length;i++){
		pcell=web_dual[list[i]];
		pcell->area=area(pcell);
		nbneighb=pcell->nb_vertices;

		//add on pressure terms
		p=pressure(pcell);
		sig.c[0][0]-=p*pcell->area;
		sig.c[1][1]-=p*pcell->area;

		div+=pcell->area;//compute divisor terms

		for(j=0;j<nbneighb;j++){
			pbond = pcell->bondlist[j];
			t = tension(pbond);
			dx=distx(pcell->vertexlist[(j+1)%nbneighb],pcell->vertexlist[j]);
			dy=disty(pcell->vertexlist[(j+1)%nbneighb],pcell->vertexlist[j]);
			d = sqrt(dx*dx+dy*dy);
			sig.c[0][0]+=t*dx*dx/d;
			sig.c[1][0]+=t*dy*dx/d;
			sig.c[0][1]+=t*dx*dy/d;
			sig.c[1][1]+=t*dy*dy/d;
		}
	}
	sig.c[0][0]/=div;sig.c[0][1]/=div;sig.c[1][0]/=div;sig.c[1][1]/=div;
	return sig;
}

matrix33 stress_cells_3d(int *list,int length){

	int i,j,nbneighb;
	Dual *pcell;
	Bond *pbond;
	double div,p,t,dx,dy,dz,d;
	matrix33 sig;

	div=0.;
	sig.c[0][0]=0.;sig.c[0][1]=0.;sig.c[1][0]=0.;sig.c[1][1]=0.;

	for(i=0;i<length;i++){
		pcell=web_dual[list[i]];
		pcell->area=area(pcell);
		nbneighb=pcell->nb_vertices;

		//add on pressure terms
		p=pressure(pcell);
		sig.c[0][0]-=p*pcell->area;
		sig.c[1][1]-=p*pcell->area;
		sig.c[2][2]-=p*pcell->area;

		div+=pcell->area;//compute divisor terms

		for(j=0;j<nbneighb;j++){
			pbond = pcell->bondlist[j];
			t = tension(pbond);
			dx=distx(pcell->vertexlist[(j+1)%nbneighb],pcell->vertexlist[j]);
			dy=disty(pcell->vertexlist[(j+1)%nbneighb],pcell->vertexlist[j]);
			dz=distz(pcell->vertexlist[(j+1)%nbneighb],pcell->vertexlist[j]);
			d = sqrt(dx*dx+dy*dy+dz*dz);
			sig.c[0][0]+=t*dx*dx/d;
			sig.c[1][0]+=t*dy*dx/d;
			sig.c[2][0]+=t*dz*dx/d;
			sig.c[0][1]+=t*dx*dy/d;
			sig.c[1][1]+=t*dy*dy/d;
			sig.c[2][1]+=t*dz*dy/d;
			sig.c[0][2]+=t*dx*dz/d;
			sig.c[1][2]+=t*dy*dz/d;
			sig.c[2][2]+=t*dz*dz/d;
			
		}
	}
	sig.c[0][0]/=div;sig.c[0][1]/=div;sig.c[1][0]/=div;sig.c[2][0]/=div;sig.c[1][1]/=div;
	sig.c[2][1]/=div;sig.c[0][2]/=div;sig.c[1][2]/=div;sig.c[2][2]/=div;

	return sig;
}
