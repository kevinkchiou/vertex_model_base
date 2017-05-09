#ifndef SRC_MEASUREMENTS_MEASUREMENTS_H
#define SRC_MEASUREMENTS_MEASUREMENTS_H 1


#define NB_COLUMNS 200
#define RAD_INC 0.05

/* measure.c */
double dist(Node *pvertex1, Node *pvertex2);

//different types of overall area calculations.  both call area_calc() below
double area(Dual * pcell);
double part_area(Dual *pcell,int max,int ivert0cell0,double x,double y,double z);
//calculates the area of the triangle created by pv1 and pv2 about centroid
vector3 area_calc(Node *pv1, Node *pv2, vector3 centroid);

//these calculate distances between different objects. useful for toroidal case
double distx(Node *pvertex1, Node *pvertex2);
double disty(Node *pvertex1, Node *pvertex2);
double distz(Node *pvertex1, Node *pvertex2);
double centdistx(vector3 centroid, Node *pvertex);
double centdisty(vector3 centroid, Node *pvertex);
double centdistz(vector3 centroid, Node *pvertex);
double ptptdistx(double x1, double x2);
double ptptdisty(double y1, double y2);
double ptptdistz(double z1, double z2);
double fixxloc(double xposition,double xbound);
double fixyloc(double yposition,double ybound);

vector3 avgpos(Node *pvert1, Node *pvert2);

vector3 centroid(Dual *pcell);

double pressure(Dual *pcell);
double perimeter(Dual *pcell);
double tension(Bond *pbond);
vector3 center_lattice();

matrix22 stress_cells_2d(int *celllist,int listlength);
matrix33 stress_cells_3d(int *celllist,int listlength);
    
#endif /* SRC_MEASUREMENTS_MEASUREMENTS_H */
