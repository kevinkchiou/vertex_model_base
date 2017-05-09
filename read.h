#ifndef SRC_DISKACCESS_READ_H
#define SRC_DISKACCESS_READ_H 1

/* Open a lattice file (ask for the file in the terminal) */
int readweb();

void generatedata(FILE *file1, FILE *file2); //basically replaces coord_gen()
void clear(); //clears the extraneous data given by coord_gen
void datadump(); //a test function to see what the hell is going wrong
void infoinput(char *vertfilename, char *cellfilename,int torus);//takes in filenames for input
#endif /* SRC_DISKACCESS_READ_H */
