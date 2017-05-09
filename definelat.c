
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>

#include "const.h" /* Main file where the parameters are defined */

#include "definelat.h"
#include "lattice.h" /* Declaration of the structures of the lattice and declaration of the functions defined in this file */
#include "locerror.h" /* Error function declaration */
#include "measurements.h" /* Declaration of the functions used for statistics */
#include "pngwrite.h"
#include "opti.h"
#include "compalone.h"
#include "in_development.h"

/* Time is evaluated from effdivision() function */
long int itime=0;
double ttime;

/* # of t1 processes */
int it1process=0;

int last_first_lattice; /* variable used to build the regular lattice */
int nb_cell_tot; /* Total number of cells (including the exterior cell which is composed of the vertices on the boundary) */
int nb_clone_tot;
int nb_vertex_tot; /* Total number of vertices */
int nb_cell_killed; 
int nb_bond_tot;

//these are non-toroidal parameters. if a toroidal lattice is used these should be changed in initialization
double LX=1000.;
double LY=1000.;
double LXMIN=-1000.;
double LYMIN=-1000.;
int TOROIDAL=0;

const int NUMDIM=2; //should be set to 2 or 3 spatial dimensions

char AREAMODELTYPE[64];
char PERIMMODELTYPE[64];

int flagt1;
int flagt1vertex;
int flagt1neighb;

double **forcevector;

/* Defines the division rate of clones and background :
 * 1 (or >1) is for a uniform division_rate
 * 0.5 makes it divide twice slower, etc
 * May be used in compalone.c
 */
double division_rates[2]={1, 1};
double cell_sizes[2]={1, 1};


unsigned int rand_disp=0;

unsigned int edges_dist[10]={0};

double add_factor = 0.0000;
double add_factor1 = 0.0000;
double add_factor2 = 0.0;

/* Pointer of the cell in the center of the initial lattice
*/
Dual *pcentercell;

/* Arrays of pointers to vertices or cells
*/
Node **web=NULL;
Dual **web_dual=NULL;
Bond **web_bond=NULL;

/* Structure used to decide whether a t1process is necessary
*/
struct t1struct {
   Node *pvertex;
   int ineighb;
} t1data;


/* Allocation functions
*/

/* Returns a pointer to a memory block of size dimension*(Bond*)*/
Bond **createbondvector(int dimension){
   Bond **a;
   a = (Bond**)malloc(dimension * sizeof(Bond*));
   if (!a) locerror("createbondvector", "Not enough memory");
   return a;
}

/* Reallocation (to different size)*/
Bond **reallocbondvector(Bond **vectortorealloc, int newdimension){
   Bond **a;
   a = (Bond**)realloc(vectortorealloc, newdimension*sizeof(Bond*));
   if (!a) locerror("reallocbondvector", "Not enough memory");
   return a;
}

/* Returns a pointer to a memory block of size dimension*(Node*)
*/
   Node
**createnodevector(int dimension)
{
   Node**a;
   a = (Node**)malloc(dimension * sizeof(Node*));
   if (!a) locerror("createnodevector", "Not enough memory");
   return a;
}

/* Reallocation (to a bigger size for example)
*/
   Node
**reallocnodevector(Node **vectortorealloc, int newdimension)
{
   Node**a;
   a = (Node**)realloc(vectortorealloc, newdimension * sizeof(Node*));
   if (!a) locerror("reallocnodevector", "Not enough memory");
   return a;
}

/* The same for cells
*/
   Dual
**createdualvector(int dimension)
{
   Dual **a;
   a = (Dual **)malloc(dimension * sizeof(Dual *));
   if (!a) locerror("createdualvector", "Not enough memory");
   return a;
}

   Dual
**reallocdualvector(Dual **vectortorealloc, int newdimension)
{
   Dual **a;
   a = (Dual **)realloc(vectortorealloc, newdimension * sizeof(Dual *));
   if (!a) locerror("reallocdualvector", "Not enough memory");
   return a;
}

/* Allocation of space for a vertex
*/
   Node
*allocnode()
{
   Node*a;
   a = (Node*)malloc(sizeof(Node));
   if (!a) locerror("allocnode", "Not enough memory");
   return a;
}

/* Allocation of space for a cell
*/
   Dual
*allocdual()
{
   Dual *a;
   a = (Dual *)malloc(sizeof(Dual));
   if (!a) locerror("allocdual", "Not enough memory");
   return a;
}

/* Allocation of space for a bond*/
Bond *allocbond(){
   Bond *a;
   a = (Bond *)malloc(sizeof(Bond));
   if (!a) locerror("allocbond", "Not enough memroy");
   return a;
}

/* Allacation of space for a list of vertices in a cell
*/
   Node
**allocvertexlist(int dimension)
{
   Node **a;
   a = (Node **)malloc(dimension * sizeof(Node *));
   if (!a) locerror("allocvertexlist", "Not enough memory");
   return a;
}

   Node
** reallocvertexlist(Node **listtorealloc, int newdimension)
{
   Node **a;
   a = (Node **)realloc(listtorealloc, newdimension * sizeof(Node *));
   if (!a) locerror("reallocvertexlist", "Not enough memory");
   return a;
}

/*Allocation of space for a list of neighboring cells*/
Dual **alloccelllist(int dimension){
   Dual **a;
   a = (Dual**)malloc(dimension*sizeof(Dual*));
   if (!a) locerror("alloccelllist", "Not enough memory");
   return a;
}

Dual **realloccelllist(Dual **listtorealloc, int newdimension){
   Dual **a;
   a = (Dual**)realloc(listtorealloc, newdimension*sizeof(Dual*));
   if (!a) locerror("realloccelllist", "Not enough memory");
   return a;
}

/*Allocation of space for list of boundary bonds*/
Bond **allocbondlist(int dimension){
   Bond **a;
   a = (Bond**)malloc(dimension * sizeof(Bond*));
   if (!a) locerror("allocbondlist","Not enough memory");
   return a;
}

Bond **reallocbondlist(Bond **listtorealloc, int newdimension){
   Bond **a;
   a = (Bond**)realloc(listtorealloc,newdimension*sizeof(Bond*));
   if (!a) locerror("reallocbondlist","Not enough memory");
   return a;
}
/* Free the memory used for a cell */
void freecell(Dual *pcell)
{
   free(pcell->vertexlist);pcell->vertexlist=NULL;
   free(pcell->celllist);pcell->celllist=NULL;
   free(pcell->bondlist);pcell->bondlist=NULL;
   free(pcell);
   pcell=NULL;
}
