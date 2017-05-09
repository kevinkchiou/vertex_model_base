#ifndef SRC_DISKACCESS_SAVE_H
#define SRC_DISKACCESS_SAVE_H 1

#include "lattice.h"
#include "definelat.h"

typedef struct nodeint
{
  double x;
  double y;		// position of the vertex 
  double z;

  int ineighb[3];		// # of the neighbours 
  int incell[3];		// # of the neighborhood cells 
} Nodeint;

typedef struct dualint
{
	int nb_vertices;
	struct strmarker marker;
} Dualint;


/* Save the current lattice in a file */
int saveweb(char *filename);
// Saves lattice to file that can be viewed by kevin's vis opengl program
void latticeprint(char *filename, unsigned int type);
// Saves lattice to *info files (useful informative format)
void infoprint(char *filename, unsigned int type);

FILE *openfile(char *filename);

int createFilenames(size_t n,char *svgfname,char *infofname,char *svgdir,char *infodir,char *fname);

#endif /* SRC_DISKACCESS_SAVE_H */
