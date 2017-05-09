#ifndef SRC_LATTICE_LATTICE_H
#define SRC_LATTICE_LATTICE_H 1

#include "const.h" /* Many preprocessor macros are used in this file */
#include "definelat.h"

/*********************************************************************/
/* We expect this new code to have 3d coordinates, but it is written */
/* with the expectation that the lattice structure lies primarily in */
/* 2d. So coords are 3d, but cells are coded considering 2d structure*/
/*********************************************************************/


/* See document file */
int vertex_from_coord(int row, int column);

void init_lattice();
void init_lattice_torus();

void kill_lattice();

void shrink_datafirst();
void shrink_data();

void cell0gen(Node *pfirstvertex);

int indexcell(Node *pvertex, Dual *pcell);
int indexvertex(Node *pvertex, Dual *pcell);

int indexvertexinweb(Node *pvertex);
int indexcellinwebdual(Dual *pcell);
int indexbondinwebbond(Bond *pbond);

void insvertex(Dual *pcell, Node *pnewvert, Node *pnextvert);
void delvertex(Dual *pcell, Node *pdelvert);

int kill_cell_p(Dual *pcell);
int kill_cell_area(Dual *pcell);
int kill_cell(Dual *pcell, int ivertex);

void division_cell(Dual *pcell, int ivertex);
void division_eq(Node *pvertex, int icell);
void division(Node *pvertex, int icell);

// These create the extended data structures for *info files
void countbonds();//counts the # of bonds
void initbonds();//initializes the bond structure
void nodeneighb();//finds neighbors to cells
void bondneighb();//finds neighbors to bonds
void neighbcell();//finds neighboring cells to cells
void createcellbondneighbs();//calls all previous functions
void updatecellbondneighbs();//same as create, but with a cleanup call
void cleanup();//cleans up the new data structures
void updatedata();//legacy version of updatecellbondneighbs();

void toruspositioning();//corrects all vertex positions on torus

void t1transform(Node *pvertex, int neighbor);
void t1flagcheck(int vertexindex, int neighbindex);//prevent multiple t1transforms


#endif /* SRC_LATTICE_LATTICE_H */
