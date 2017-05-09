
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>

#include "const.h" /* Main file where the parameters are defined */

#include "generatelat.h"
#include "definelat.h"
#include "lattice.h" /* Declaration of the structures of the lattice and declaration of the functions defined in this file */
#include "locerror.h" /* Error function declaration */
#include "measurements.h" /* Declaration of the functions used for statistics */
#include "pngwrite.h"
#include "opti.h"
#include "compalone.h"
#include "in_development.h"

/* Generation of the first lattice (some pointers to cells remains void and shrink_datafirst() must be applied)
*/
void coord_gen_torus2(){ //creates torus with x-axis twice as long as other version
    int n,i,k,l,m,kk;
    double sq3=sqrt(3.),xc,yc,cox,coy,tol=0.001,LXLOC,LYLOC;
    Dual *pcell;

    LX=40.;
    LY=LX*sqrt(3.)*0.25;
    LXMIN=0.;
    LYMIN=0.; //setting the toroidal boundary
    TOROIDAL=1;


    nb_vertex_tot=768;
    nb_cell_tot=384;

    //nb_vertex_tot=192;
    //nb_cell_tot=96;

    //nb_vertex_tot=48;
    //nb_cell_tot=24;

    LXLOC=sqrt(6*nb_cell_tot);
    LYLOC=LXLOC*0.25*sqrt(3.);


    web      = createnodevector(nb_vertex_tot);
    web_dual = createdualvector(nb_cell_tot);

    for (n=0;n<nb_vertex_tot;n++){
        web[n]=allocnode();
    }
    for (n=0;n<nb_cell_tot;n++){
        web_dual[n]=allocdual();
    }

    //Everything is counter clock wise
    //Positions of vertices:

    i=0;
    for (k=0;k<=(int)((LXLOC-0.5)/3.);k++){
        for (l=0;l<=(int)(LYLOC/sq3-0.5);l++){
            web[i]->x=0.5+3*k;
            web[i]->y=0.5*sq3*(2*l+1);
            i++;
        }
    }
    for (k=0;k<=(int)((LXLOC-2)/3.);k++){
        for (l=0;l<(int)(LYLOC/sq3);l++){
            web[i]->x=2.+3*k;
            web[i]->y=sq3*l;
            i++;
        }
    }

    for (k=0;k<=(int)((LXLOC-1)/3.);k++){
        for (l=0;l<(int)(LYLOC/sq3);l++){
            web[i]->x=1+3*k;
            web[i]->y=sq3*l;
            i++;
        }
    }
    for (k=0;k<=(int)((LXLOC-2.5)/3.);k++){
        for (l=0;l<=(int)(LYLOC/sq3-0.5);l++){
            web[i]->x=2.5+3*k;
            web[i]->y=0.5*sq3*(2*l+1);
            i++;
        }
    }

    //Vertex neighbors  

    for (i=0;i<nb_vertex_tot;i++){

        cox = fixxloc(web[i]->x-1,LXLOC);
        coy = fixyloc(web[i]->y,LYLOC);

        //Neighbor 0
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol ) && (fabs(web[k]->y-coy)<tol )){
                web[i]->pneighb[0]=web[k];
                break;
            }
        }

        cox = fixxloc(web[i]->x+1,LXLOC);
        coy = fixyloc(web[i]->y,LYLOC);
        //Neighbor 0
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol) && (fabs(web[k]->y-coy)<tol)){
                web[i]->pneighb[0]=web[k];
                break;
            }
        }
        //Neighbor 1
        cox = fixxloc(web[i]->x+0.5,LXLOC);
        coy = fixyloc(web[i]->y-0.5*sq3,LYLOC);
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol)){
                web[i]->pneighb[1]=web[k];
                break;
            }
        }
        //Neighbor 1
        cox = fixxloc(web[i]->x-0.5,LXLOC);
        coy = fixyloc(web[i]->y+0.5*sq3,LYLOC);
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol)){
                web[i]->pneighb[1]=web[k];
                break;
            }
        }
        //Neighbor 2
        cox = fixxloc(web[i]->x+0.5,LXLOC);
        coy = fixyloc(web[i]->y+0.5*sq3,LYLOC);
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol)){
                web[i]->pneighb[2]=web[k];
                break;
            }
        }
        //Neighbor 2
        cox = fixxloc(web[i]->x-0.5,LXLOC);
        coy = fixyloc(web[i]->y-0.5*sq3,LYLOC);
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol)){
                web[i]->pneighb[2]=web[k];
                break;
            }
        }
    }

    //   for (k=0;k<nb_vertex_tot;k++){ 
    //     for (i=0;i<3;i++){ 
    //       printf("%i %i %i\n",k,i,indexvertexinweb(web[k]->pneighb[i])); 
    //     }
    //   }
    //   printf("done with the neighbors!\n");fflush(stdout);

    //Cells
    i=0;

    for (kk=0;kk<(int)(LXLOC/3.);kk++){
        for (l=0;l<(int)(LYLOC/sq3);l++){
            web_dual[i]->nb_vertices = 6;
            web_dual[i]->vertexlist=allocvertexlist(web_dual[i]->nb_vertices);
            xc=3.*kk;
            yc=sq3*l;
            //Neighbor 0
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-1.,LXLOC);
                coy=fixyloc(yc,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[0]=web[k];
                    break;
                }
            }
            //Neighbor 1
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-0.5,LXLOC);
                coy=fixyloc(yc-sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[1]=web[k];
                    break;
                }
            }
            //Neighbor 2
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+0.5,LXLOC);
                coy=fixyloc(yc-sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[2]=web[k];
                    break;
                }
            }
            //Neighbor 3
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+1,LXLOC);
                coy=fixyloc(yc,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[3]=web[k];
                    break;
                }
            }
            //Neighbor 4
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+0.5,LXLOC);
                coy=fixyloc(yc+sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[4]=web[k];
                    break;
                }
            }
            //Neighbor 5
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-0.5,LXLOC);
                coy=fixyloc(yc+sq3/2,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[5]=web[k];
                    break;
                }
            }
            i++;
        }
    }


    for (kk=0;kk<=(int)((LXLOC-1.5)/3.);kk++){
        for (l=0;l<=(int)(LYLOC/sq3-0.5);l++){
            web_dual[i]->nb_vertices = 6;
            web_dual[i]->vertexlist=allocvertexlist(web_dual[i]->nb_vertices);
            xc=1.5+3*kk;
            yc=0.5*sq3*(2*l+1);
            //Neighbor 0
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-1.,LXLOC);
                coy=fixyloc(yc,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[0]=web[k];
                    break;
                }
            }
            //Neighbor 1
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-0.5,LXLOC);
                coy=fixyloc(yc-sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[1]=web[k];
                    break;
                }
            }
            //Neighbor 2
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+0.5,LXLOC);
                coy=fixyloc(yc-sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[2]=web[k];
                    break;
                }
            }
            //Neighbor 3
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+1,LXLOC);
                coy=fixyloc(yc,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[3]=web[k];
                    break;
                }
            }
            //Neighbor 4
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+0.5,LXLOC);
                coy=fixyloc(yc+sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[4]=web[k];
                    break;
                }
            }
            //Neighbor 5
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-0.5,LXLOC);
                coy=fixyloc(yc+sq3/2,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[5]=web[k];
                    break;
                }
            }
            i++;

        }
    }


    //Neighbor cells


    for (i=0;i<nb_vertex_tot;i++){
        for (l=0;l<3;l++){
            for (k=0;k<nb_cell_tot;k++){
                for (m=0;m<web_dual[k]->nb_vertices;m++){
                    if ((web_dual[k]->vertexlist[m]==web[i]->pneighb[(l+1)%3]) && (web_dual[k]->vertexlist[(m+1)%web_dual[k]->nb_vertices]==web[i]) && (web_dual[k]->vertexlist[(m+2)%web_dual[k]->nb_vertices]==web[i]->pneighb[l])) web[i]->pncell[l]=web_dual[k];
                }
            }
        }
    }

    for (i=0;i<nb_vertex_tot;i++){
        web[i]->x*=LX/LXLOC;
        web[i]->y*=LY/LYLOC;
    }
}

void coord_gen_torus(){
    int n,i,k,l,m,kk;
    double sq3=sqrt(3.),xc,yc,cox,coy,tol=0.001,LXLOC,LYLOC;
    Dual *pcell;

    LX=20.;
    LY=LX*sqrt(3.)*0.5;
    LXMIN=0.;
    LYMIN=0.; //setting the toroidal boundary
    TOROIDAL=1;


    nb_vertex_tot=384;
    nb_cell_tot=192;

    //nb_vertex_tot=96;
    //nb_cell_tot=48;

    //nb_vertex_tot=24;
    //nb_cell_tot=12;

    LXLOC=sqrt(3*nb_cell_tot);
    LYLOC=LXLOC*0.5*sqrt(3.);


    web      = createnodevector(nb_vertex_tot);
    web_dual = createdualvector(nb_cell_tot);

    for (n=0;n<nb_vertex_tot;n++){
        web[n]=allocnode();
    }
    for (n=0;n<nb_cell_tot;n++){
        web_dual[n]=allocdual();
    }

    //Everything is counter clock wise
    //Positions of vertices:

    i=0;
    for (k=0;k<=(int)((LXLOC-0.5)/3.);k++){
        for (l=0;l<=(int)(LYLOC/sq3-0.5);l++){
            web[i]->x=0.5+3*k;
            web[i]->y=0.5*sq3*(2*l+1);
            i++;
        }
    }
    for (k=0;k<=(int)((LXLOC-2)/3.);k++){
        for (l=0;l<(int)(LYLOC/sq3);l++){
            web[i]->x=2.+3*k;
            web[i]->y=sq3*l;
            i++;
        }
    }

    for (k=0;k<=(int)((LXLOC-1)/3.);k++){
        for (l=0;l<(int)(LYLOC/sq3);l++){
            web[i]->x=1+3*k;
            web[i]->y=sq3*l;
            i++;
        }
    }
    for (k=0;k<=(int)((LXLOC-2.5)/3.);k++){
        for (l=0;l<=(int)(LYLOC/sq3-0.5);l++){
            web[i]->x=2.5+3*k;
            web[i]->y=0.5*sq3*(2*l+1);
            i++;
        }
    }

    //Vertex neighbors  

    for (i=0;i<nb_vertex_tot;i++){

        cox = fixxloc(web[i]->x-1,LXLOC);
        coy = fixyloc(web[i]->y,LYLOC);

        //Neighbor 0
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol ) && (fabs(web[k]->y-coy)<tol )){
                web[i]->pneighb[0]=web[k];
                break;
            }
        }

        cox = fixxloc(web[i]->x+1,LXLOC);
        coy = fixyloc(web[i]->y,LYLOC);
        //Neighbor 0
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol) && (fabs(web[k]->y-coy)<tol)){
                web[i]->pneighb[0]=web[k];
                break;
            }
        }
        //Neighbor 1
        cox = fixxloc(web[i]->x+0.5,LXLOC);
        coy = fixyloc(web[i]->y-0.5*sq3,LYLOC);
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol)){
                web[i]->pneighb[1]=web[k];
                break;
            }
        }
        //Neighbor 1
        cox = fixxloc(web[i]->x-0.5,LXLOC);
        coy = fixyloc(web[i]->y+0.5*sq3,LYLOC);
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol)){
                web[i]->pneighb[1]=web[k];
                break;
            }
        }
        //Neighbor 2
        cox = fixxloc(web[i]->x+0.5,LXLOC);
        coy = fixyloc(web[i]->y+0.5*sq3,LYLOC);
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol)){
                web[i]->pneighb[2]=web[k];
                break;
            }
        }
        //Neighbor 2
        cox = fixxloc(web[i]->x-0.5,LXLOC);
        coy = fixyloc(web[i]->y-0.5*sq3,LYLOC);
        for (k=0;k<nb_vertex_tot;k++){
            if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol)){
                web[i]->pneighb[2]=web[k];
                break;
            }
        }
    }

    //   for (k=0;k<nb_vertex_tot;k++){ 
    //     for (i=0;i<3;i++){ 
    //       printf("%i %i %i\n",k,i,indexvertexinweb(web[k]->pneighb[i])); 
    //     }
    //   }
    //   printf("done with the neighbors!\n");fflush(stdout);

    //Cells
    i=0;

    for (kk=0;kk<(int)(LXLOC/3.);kk++){
        for (l=0;l<(int)(LYLOC/sq3);l++){
            web_dual[i]->nb_vertices = 6;
            web_dual[i]->vertexlist=allocvertexlist(web_dual[i]->nb_vertices);
            xc=3.*kk;
            yc=sq3*l;
            //Neighbor 0
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-1.,LXLOC);
                coy=fixyloc(yc,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[0]=web[k];
                    break;
                }
            }
            //Neighbor 1
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-0.5,LXLOC);
                coy=fixyloc(yc-sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[1]=web[k];
                    break;
                }
            }
            //Neighbor 2
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+0.5,LXLOC);
                coy=fixyloc(yc-sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[2]=web[k];
                    break;
                }
            }
            //Neighbor 3
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+1,LXLOC);
                coy=fixyloc(yc,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[3]=web[k];
                    break;
                }
            }
            //Neighbor 4
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+0.5,LXLOC);
                coy=fixyloc(yc+sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[4]=web[k];
                    break;
                }
            }
            //Neighbor 5
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-0.5,LXLOC);
                coy=fixyloc(yc+sq3/2,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[5]=web[k];
                    break;
                }
            }
            i++;
        }
    }


    for (kk=0;kk<=(int)((LXLOC-1.5)/3.);kk++){
        for (l=0;l<=(int)(LYLOC/sq3-0.5);l++){
            web_dual[i]->nb_vertices = 6;
            web_dual[i]->vertexlist=allocvertexlist(web_dual[i]->nb_vertices);
            xc=1.5+3*kk;
            yc=0.5*sq3*(2*l+1);
            //Neighbor 0
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-1.,LXLOC);
                coy=fixyloc(yc,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[0]=web[k];
                    break;
                }
            }
            //Neighbor 1
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-0.5,LXLOC);
                coy=fixyloc(yc-sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[1]=web[k];
                    break;
                }
            }
            //Neighbor 2
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+0.5,LXLOC);
                coy=fixyloc(yc-sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[2]=web[k];
                    break;
                }
            }
            //Neighbor 3
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+1,LXLOC);
                coy=fixyloc(yc,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[3]=web[k];
                    break;
                }
            }
            //Neighbor 4
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc+0.5,LXLOC);
                coy=fixyloc(yc+sq3/2.,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[4]=web[k];
                    break;
                }
            }
            //Neighbor 5
            for (k=0;k<nb_vertex_tot;k++){
                cox=fixxloc(xc-0.5,LXLOC);
                coy=fixyloc(yc+sq3/2,LYLOC);
                if ((fabs(web[k]->x-cox)<tol || fabs(web[k]->x+LXLOC-cox)<tol || fabs(web[k]->x-LXLOC-cox)<tol) && (fabs(web[k]->y-coy)<tol || fabs(web[k]->y+LYLOC-coy)<tol || fabs(web[k]->y-LYLOC-coy)<tol)){
                    web_dual[i]->vertexlist[5]=web[k];
                    break;
                }
            }
            i++;

        }
    }


    //Neighbor cells


    for (i=0;i<nb_vertex_tot;i++){
        for (l=0;l<3;l++){
            for (k=0;k<nb_cell_tot;k++){
                for (m=0;m<web_dual[k]->nb_vertices;m++){
                    if ((web_dual[k]->vertexlist[m]==web[i]->pneighb[(l+1)%3]) && (web_dual[k]->vertexlist[(m+1)%web_dual[k]->nb_vertices]==web[i]) && (web_dual[k]->vertexlist[(m+2)%web_dual[k]->nb_vertices]==web[i]->pneighb[l])) web[i]->pncell[l]=web_dual[k];
                }
            }
        }
    }

    for (i=0;i<nb_vertex_tot;i++){
        web[i]->x*=LX/LXLOC;
        web[i]->y*=LY/LYLOC;
    }
}

void coord_gen()
{
   int row, column, column_lim, vertex, cell;
   Node *pvertex;
   Dual *pcell;

   last_first_lattice=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1);

   web = createnodevector(last_first_lattice);

   web_dual=createdualvector(last_first_lattice);

   web_dual[0]=allocdual();
   web_dual[0] -> area =0;

   /* Allocation of memory for the vertices
    * The code is not optimized at all but it is only used for initialisation and it is more readable
    */

   for (row=2 ; row<2*LENGTH_FIRST_LATTICE-1; row++){

      if (row%2){
         column_lim=LENGTH_FIRST_LATTICE + 1;
         for (column=1 ; column<column_lim; column++){

            vertex=vertex_from_coord(row, column);
            web[vertex]=allocnode();
         }
      }

      if ((row+1)%2){
         column_lim=LENGTH_FIRST_LATTICE - 1;
         for (column=1 ; column<column_lim; column++){

            vertex=vertex_from_coord(row, column);
            web[vertex]=allocnode();
         }

         for (column=column_lim ; column < column_lim + 2; column++){
            vertex=vertex_from_coord(row, column);
            web[vertex]=NULL;
         }
      }
   }

   web[0]=NULL;
   web[1]=NULL;
   web[LENGTH_FIRST_LATTICE]=NULL;

   for (column=2 ; column < LENGTH_FIRST_LATTICE; column++){
      vertex=column;
      web[vertex]=allocnode();
   }

   web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+1]=NULL;

   for (column=2 ; column < LENGTH_FIRST_LATTICE; column++){
      vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+column;
      web[vertex]=allocnode();
   }

   /* Allocation of memory for the cells */

   for (row=1 ; row<2*LENGTH_FIRST_LATTICE-1; row++){

      if (row%2){
         column=1;
         vertex=vertex_from_coord(row, column);
         web_dual[vertex]=NULL;

         column_lim=LENGTH_FIRST_LATTICE + 1;
         web_dual[LENGTH_FIRST_LATTICE*(row-1)+1]=NULL;
         for (column=2 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            if (column%2){
               web_dual[vertex]=allocdual();
            }
            else {
               web_dual[vertex]=NULL;
            }
         }
      }

      if ((row+1)%2){
         column_lim=LENGTH_FIRST_LATTICE+1;
         for (column=1 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);

            if (column%2){
               web_dual[vertex]=allocdual();
            }
            else {
               web_dual[vertex]=NULL;
            }
         }
      }
   }
   row=2*LENGTH_FIRST_LATTICE-1;
   column=1;
   vertex=vertex_from_coord(row, column);
   web_dual[vertex]=NULL;

   column_lim=LENGTH_FIRST_LATTICE;
   for (column=2 ; column<column_lim; column++){
      vertex=vertex_from_coord(row, column);
      if (column%2){
         web_dual[vertex]=allocdual();
      }
      else {
         web_dual[vertex]=NULL;
      }
   }


   web_dual[0]=allocdual();

   /* Effective construction of the lattice */

   for (row=2 ; row < 2*LENGTH_FIRST_LATTICE -1 ; row++){

      if (row%2){

         column_lim=LENGTH_FIRST_LATTICE;

         for (column=2 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];
            if (pvertex){

               pvertex -> y=(row-1)*sqrt(3)/2;

               if (column%2){
                  pvertex -> x=-1.5+1.5*column;

                  pvertex -> pneighb[0]=web[vertex+1];
                  pvertex -> pneighb[1]=web[vertex+LENGTH_FIRST_LATTICE-1];
                  pvertex -> pneighb[2]=web[vertex-LENGTH_FIRST_LATTICE-1];

                  pvertex -> pncell[0]=web_dual[vertex+LENGTH_FIRST_LATTICE];	
                  pvertex -> pncell[1]=web_dual[vertex];
                  pvertex -> pncell[2]=web_dual[vertex-LENGTH_FIRST_LATTICE];
               }

               if ((column+1)%2){
                  pvertex -> x=-2+1.5*column;

                  pvertex -> pneighb[0]=web[vertex-1];
                  pvertex -> pneighb[1]=web[vertex-LENGTH_FIRST_LATTICE-1];
                  pvertex -> pneighb[2]=web[vertex+LENGTH_FIRST_LATTICE-1];

                  pvertex -> pncell[0]=web_dual[vertex-LENGTH_FIRST_LATTICE-1];	
                  pvertex -> pncell[1]=web_dual[vertex+1];
                  pvertex -> pncell[2]=web_dual[vertex+LENGTH_FIRST_LATTICE-1];
               }

            }
            else {
               locerror("coord_gen", "The vertex does not exist");
            }
         }
      }


      if ((row+1)%2){

         column_lim=LENGTH_FIRST_LATTICE-1;

         for (column=1 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];
            if (pvertex){

               pvertex -> y=(row-1)*sqrt(3)/2;

               if (column%2){
                  pvertex -> x=1.5*column;

                  pvertex -> pneighb[0]=web[vertex+1];
                  pvertex -> pneighb[1]=web[vertex+LENGTH_FIRST_LATTICE+1];
                  pvertex -> pneighb[2]=web[vertex-LENGTH_FIRST_LATTICE+1];

                  pvertex -> pncell[0]=web_dual[vertex+LENGTH_FIRST_LATTICE+2];	
                  pvertex -> pncell[1]=web_dual[vertex];
                  pvertex -> pncell[2]=web_dual[vertex-LENGTH_FIRST_LATTICE+2];
               }

               if ((column+1)%2){
                  pvertex -> x=-0.5+1.5*column;

                  pvertex -> pneighb[0]=web[vertex-1];
                  pvertex -> pneighb[1]=web[vertex-LENGTH_FIRST_LATTICE+1];
                  pvertex -> pneighb[2]=web[vertex+LENGTH_FIRST_LATTICE+1];

                  pvertex -> pncell[0]=web_dual[vertex-LENGTH_FIRST_LATTICE+1];	
                  pvertex -> pncell[1]=web_dual[vertex+1];
                  pvertex -> pncell[2]=web_dual[vertex+LENGTH_FIRST_LATTICE+1];
               }

            }
            else {
               locerror("coord_gen", "The vertex does not exist");
            }

         }
      }
   }

   column_lim=LENGTH_FIRST_LATTICE-1;
   for (column=3 ; column<column_lim; column++){
      if (web[column] && web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+column]){
         if (column%2){
            row=1;
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];
            pvertex -> x = -1.5 + 1.5*column;
            pvertex -> y = 0;

            pvertex -> pneighb[0]=web[vertex+1];
            pvertex -> pneighb[1]=web[vertex+LENGTH_FIRST_LATTICE-1];
            pvertex -> pneighb[2]=web[vertex-1];

            pvertex -> pncell[0]=web_dual[vertex+LENGTH_FIRST_LATTICE];	
            pvertex -> pncell[1]=web_dual[vertex];
            pvertex -> pncell[2]=web_dual[0];


            row=2*LENGTH_FIRST_LATTICE-1;
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];

            pvertex -> x = -1.5 + 1.5*column;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=web[vertex+1];
            pvertex -> pneighb[1]=web[vertex-1];
            pvertex -> pneighb[2]=web[vertex-LENGTH_FIRST_LATTICE-1];

            pvertex -> pncell[0]=web_dual[0];	
            pvertex -> pncell[1]=web_dual[vertex];
            pvertex -> pncell[2]=web_dual[vertex-LENGTH_FIRST_LATTICE];
         }

         if ((column+1)%2){
            row=1;
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];

            pvertex -> x = -2 + 1.5*column;
            pvertex -> y = 0;

            pvertex -> pneighb[0]=web[vertex-1];
            pvertex -> pneighb[1]=web[vertex+1];
            pvertex -> pneighb[2]=web[vertex+LENGTH_FIRST_LATTICE-1];

            pvertex -> pncell[0]=web_dual[0];	
            pvertex -> pncell[1]=web_dual[vertex+1];
            pvertex -> pncell[2]=web_dual[vertex + LENGTH_FIRST_LATTICE-1];


            row=2*LENGTH_FIRST_LATTICE-1;
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];


            pvertex -> x = -2 + 1.5*column;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=web[vertex-1];
            pvertex -> pneighb[1]=web[vertex-LENGTH_FIRST_LATTICE-1];
            pvertex -> pneighb[2]=web[vertex+1];

            pvertex -> pncell[0]=web_dual[vertex-LENGTH_FIRST_LATTICE-1];
            pvertex -> pncell[1]=web_dual[vertex+1];
            pvertex -> pncell[2]=web_dual[0];
         }
      }
      else {
         locerror("coord_gen", "The vertex does not exist");
      }
   }

   for (row=5; row< 2 * LENGTH_FIRST_LATTICE - 4; row++){

      if (row%2){
         column=1;
         vertex=vertex_from_coord(row, column);
         pvertex=web[vertex];

         if (pvertex){

            pvertex -> x = 0;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=web[vertex+1];
            pvertex -> pneighb[1]=web[vertex+2*LENGTH_FIRST_LATTICE];
            pvertex -> pneighb[2]=web[vertex-2*LENGTH_FIRST_LATTICE];

            pvertex -> pncell[0]=web_dual[vertex+LENGTH_FIRST_LATTICE];
            pvertex -> pncell[1]=web_dual[0];
            pvertex -> pncell[2]=web_dual[vertex-LENGTH_FIRST_LATTICE];
         }
         else {
            locerror("coord_gen", "The vertex does not exist");
         }

         column=LENGTH_FIRST_LATTICE;
         vertex=vertex_from_coord(row, column);
         pvertex=web[vertex];

         if (pvertex){

            pvertex -> x = -2 + 1.5*LENGTH_FIRST_LATTICE;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=web[vertex-1];
            pvertex -> pneighb[1]=web[vertex-2*LENGTH_FIRST_LATTICE];
            pvertex -> pneighb[2]=web[vertex+2*LENGTH_FIRST_LATTICE];

            pvertex -> pncell[0]=web_dual[vertex-LENGTH_FIRST_LATTICE-1];
            pvertex -> pncell[1]=web_dual[0];
            pvertex -> pncell[2]=web_dual[vertex+LENGTH_FIRST_LATTICE-1];
         }
         else {
            locerror("coord_gen", "The vertex does not exist");
         }
      }
   }

   /* Corner generation */
   vertex=2;

   pvertex=web[vertex];

   pvertex -> x = 1;
   pvertex -> y = 0;

   pvertex -> pneighb[0]=web[2*LENGTH_FIRST_LATTICE+1];
   pvertex -> pneighb[1]=web[3];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE+1];

   pvertex -> pncell[0]=web_dual[0];
   pvertex -> pncell[1]=web_dual[3];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE + 1];


   vertex=2*LENGTH_FIRST_LATTICE+1;
   pvertex=web[vertex];

   pvertex -> x = 0;
   pvertex -> y = sqrt(3);

   pvertex -> pneighb[0]=web[2*LENGTH_FIRST_LATTICE+2];
   pvertex -> pneighb[1]=web[4*LENGTH_FIRST_LATTICE+1];
   pvertex -> pneighb[2]=web[2];

   pvertex -> pncell[0]=web_dual[3*LENGTH_FIRST_LATTICE+1];
   pvertex -> pncell[1]=web_dual[0];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE+1];


   vertex=LENGTH_FIRST_LATTICE-1;
   pvertex=web[vertex];

   pvertex -> x = -3 + 1.5*LENGTH_FIRST_LATTICE;
   pvertex -> y = 0;

   pvertex -> pneighb[0]=web[3*LENGTH_FIRST_LATTICE];
   pvertex -> pneighb[1]=web[2*LENGTH_FIRST_LATTICE-2];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE-2];

   pvertex -> pncell[0]=web_dual[2*LENGTH_FIRST_LATTICE-1];
   pvertex -> pncell[1]=web_dual[LENGTH_FIRST_LATTICE-1];
   pvertex -> pncell[2]=web_dual[0];


   vertex=3*LENGTH_FIRST_LATTICE;
   pvertex=web[vertex];

   pvertex -> x = -2 + 1.5*LENGTH_FIRST_LATTICE;
   pvertex -> y = sqrt(3);

   pvertex -> pneighb[0]=web[3*LENGTH_FIRST_LATTICE-1];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE-1];
   pvertex -> pneighb[2]=web[5*LENGTH_FIRST_LATTICE];

   pvertex -> pncell[0]=web_dual[2*LENGTH_FIRST_LATTICE-1];
   pvertex -> pncell[1]=web_dual[0];
   pvertex -> pncell[2]=web_dual[4*LENGTH_FIRST_LATTICE-1];


   vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+1;
   pvertex=web[vertex];

   pvertex -> x = 0;
   pvertex -> y = (2*LENGTH_FIRST_LATTICE-4)*sqrt(3)/2;

   pvertex -> pneighb[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+2];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+2];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-6)+1];

   pvertex -> pncell[0]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)+1];
   pvertex -> pncell[1]=web_dual[0];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-5)+1];


   vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+2;
   pvertex=web[vertex];

   pvertex -> x = 1;
   pvertex -> y = (2*LENGTH_FIRST_LATTICE-2)*sqrt(3)/2;

   pvertex -> pneighb[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+1];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)+1];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+3];

   pvertex -> pncell[0]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)+1];
   pvertex -> pncell[1]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+3];
   pvertex -> pncell[2]=web_dual[0];


   vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3);
   pvertex=web[vertex];

   pvertex -> x = -2 + 1.5*LENGTH_FIRST_LATTICE;
   pvertex -> y = (2*LENGTH_FIRST_LATTICE-4)*sqrt(3)/2;

   pvertex -> pneighb[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)-1];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-5)];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-1];

   pvertex -> pncell[0]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)-1];
   pvertex -> pncell[1]=web_dual[0];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-1];


   vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-1;
   pvertex=web[vertex];

   pvertex -> x = -3 + 1.5*LENGTH_FIRST_LATTICE;
   pvertex -> y = (2*LENGTH_FIRST_LATTICE-2)*sqrt(3)/2;

   pvertex -> pneighb[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-2];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-2];

   pvertex -> pncell[0]=web_dual[0];
   pvertex -> pncell[1]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-1];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-1];

   /* List of vertices in cells */

   /* Hexagons */

   for (row=2 ; row<2*LENGTH_FIRST_LATTICE-1; row++){

      if (row%2){
         column_lim=LENGTH_FIRST_LATTICE;
         for (column=3 ; column<column_lim; column++){
            if (column%2){

               cell=vertex_from_coord(row, column);
               pcell=web_dual[cell];
               if (pcell){
                  pcell -> nb_vertices = 6;
                  pcell -> vertexlist = allocvertexlist(6);
                  pcell -> vertexlist[0]=web[cell];
                  pcell -> vertexlist[1]=web[cell+LENGTH_FIRST_LATTICE-1];
                  pcell -> vertexlist[2]=web[cell+LENGTH_FIRST_LATTICE-2];
                  pcell -> vertexlist[3]=web[cell-1];
                  pcell -> vertexlist[4]=web[cell-LENGTH_FIRST_LATTICE-2];
                  pcell -> vertexlist[5]=web[cell-LENGTH_FIRST_LATTICE-1];
               }
               else {
                  locerror("coord_gen", "The cell does not exist");
               }
            }

         }
      }

      if ((row+1)%2){
         column_lim=LENGTH_FIRST_LATTICE - 2;
         for (column=3 ; column<column_lim; column++){

            if (column%2){

               cell=vertex_from_coord(row, column);
               pcell=web_dual[cell];
               if (pcell){
                  pcell -> nb_vertices = 6;
                  pcell -> vertexlist = allocvertexlist(6);
                  pcell -> vertexlist[0]=web[cell];
                  pcell -> vertexlist[1]=web[cell+LENGTH_FIRST_LATTICE+1];
                  pcell -> vertexlist[2]=web[cell+LENGTH_FIRST_LATTICE];
                  pcell -> vertexlist[3]=web[cell-1];
                  pcell -> vertexlist[4]=web[cell-LENGTH_FIRST_LATTICE];
                  pcell -> vertexlist[5]=web[cell-LENGTH_FIRST_LATTICE+1];
               }
               else {
                  locerror("coord_gen", "The cell does not exist");
               }
            }
         }
      }
   }

   for (row=4 ; row<2*LENGTH_FIRST_LATTICE-3; row++){

      if ((row+1)%2){
         column=1;
         cell=vertex_from_coord(row, column);
         pcell=web_dual[cell];

         if (pcell){
            pcell -> nb_vertices = 5;
            pcell -> vertexlist = allocvertexlist(5);
            pcell -> vertexlist[0]=web[cell];
            pcell -> vertexlist[1]=web[cell+LENGTH_FIRST_LATTICE+1];
            pcell -> vertexlist[2]=web[cell+LENGTH_FIRST_LATTICE];
            pcell -> vertexlist[3]=web[cell-LENGTH_FIRST_LATTICE];
            pcell -> vertexlist[4]=web[cell-LENGTH_FIRST_LATTICE+1];
         }
         else {
            locerror("coord_gen", "The cell does not exist");
         }


         column=LENGTH_FIRST_LATTICE-2;
         cell=vertex_from_coord(row, column);
         pcell=web_dual[cell+1];

         if (pcell){
            pcell -> nb_vertices = 5;
            pcell -> vertexlist = allocvertexlist(5);
            pcell -> vertexlist[0]=web[cell];
            pcell -> vertexlist[1]=web[cell-LENGTH_FIRST_LATTICE+1];
            pcell -> vertexlist[2]=web[cell-LENGTH_FIRST_LATTICE+2];
            pcell -> vertexlist[3]=web[cell+LENGTH_FIRST_LATTICE+2];
            pcell -> vertexlist[4]=web[cell+LENGTH_FIRST_LATTICE+1];
         }
         else {
            locerror("coord_gen", "The cell does not exist");
         }

      }
   }

   column_lim=LENGTH_FIRST_LATTICE-1;
   for (column=2 ; column<column_lim; column++){
      if ((column+1)%2){
         row=1;
         cell=vertex_from_coord(row, column);
         pcell=web_dual[cell+1];
         if (pcell){
            pcell -> nb_vertices = 4;
            pcell -> vertexlist = allocvertexlist(4);
            pcell -> vertexlist[0]=web[cell];
            pcell -> vertexlist[1]=web[cell+1];
            pcell -> vertexlist[2]=web[cell+LENGTH_FIRST_LATTICE];
            pcell -> vertexlist[3]=web[cell+LENGTH_FIRST_LATTICE-1];
         }
         else {
            locerror("coord_gen", "The cell+1 does not exist");
         }


         row=2*LENGTH_FIRST_LATTICE-1;
         cell=vertex_from_coord(row, column);
         pcell=web_dual[cell+1];
         if (pcell){
            pcell -> nb_vertices = 4;
            pcell -> vertexlist = allocvertexlist(4);
            pcell -> vertexlist[0]=web[cell];
            pcell -> vertexlist[1]=web[cell-LENGTH_FIRST_LATTICE-1];
            pcell -> vertexlist[2]=web[cell-LENGTH_FIRST_LATTICE];
            pcell -> vertexlist[3]=web[cell+1];
         }
         else {
            locerror("coord_gen", "The cell+1 does not exist");
         }
      }
   }


   pcell=web_dual[LENGTH_FIRST_LATTICE+1];
   if (pcell){
      pcell -> nb_vertices = 4;
      pcell -> vertexlist = allocvertexlist(4);
      pcell -> vertexlist[0]=web[LENGTH_FIRST_LATTICE+1];
      pcell -> vertexlist[1]=web[2*LENGTH_FIRST_LATTICE+2];
      pcell -> vertexlist[2]=web[2*LENGTH_FIRST_LATTICE+1];
      pcell -> vertexlist[3]=web[2];
   }
   else {
      locerror("coord_gen", "The LENGTH_FIRST_LATTICE+1 does not exist");
   }

   pcell=web_dual[2*LENGTH_FIRST_LATTICE-1];
   if (pcell){
      pcell -> nb_vertices = 4;
      pcell -> vertexlist = allocvertexlist(4);
      pcell -> vertexlist[0]=web[3*LENGTH_FIRST_LATTICE];
      pcell -> vertexlist[1]=web[3*LENGTH_FIRST_LATTICE-1];
      pcell -> vertexlist[2]=web[2*LENGTH_FIRST_LATTICE-2];
      pcell -> vertexlist[3]=web[LENGTH_FIRST_LATTICE-1];
   }
   else {
      locerror("coord_gen", "The LENGTH_FIRST_LATTICE+1 does not exist");
   }

   cell=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)+1;
   pcell=web_dual[cell];
   if (pcell){
      pcell -> nb_vertices = 4;
      pcell -> vertexlist = allocvertexlist(4);
      pcell -> vertexlist[0]=web[cell];
      pcell -> vertexlist[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+2];
      pcell -> vertexlist[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+1];
      pcell -> vertexlist[3]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+2];
   }
   else {
      locerror("coord_gen", "The LENGTH_FIRST_LATTICE+1 does not exist");
   }

   cell=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-1;
   pcell=web_dual[cell];
   if (pcell){
      pcell -> nb_vertices = 4;
      pcell -> vertexlist = allocvertexlist(4);
      pcell -> vertexlist[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-1];
      pcell -> vertexlist[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-2];
      pcell -> vertexlist[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)-1];
      pcell -> vertexlist[3]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)];
   }
   else {
      locerror("coord_gen", "The LENGTH_FIRST_LATTICE+1 does not exist");
   }


   /* Initialisation of the outside cell */

   pcell=web_dual[0];
   pcell -> nb_vertices = 4*(LENGTH_FIRST_LATTICE-2);
   pcell -> vertexlist = allocvertexlist(4*(LENGTH_FIRST_LATTICE-2));


   for (row=3 ; row<2*LENGTH_FIRST_LATTICE-2; row++){

      if (row%2){
         column=1;
         vertex=vertex_from_coord(row, column);
         pvertex=web[vertex];
         pcell -> vertexlist[(row-1)/2-1]=pvertex;

         column=LENGTH_FIRST_LATTICE;
         vertex=vertex_from_coord(row, column);
         pvertex=web[vertex];
         pcell -> vertexlist[3*LENGTH_FIRST_LATTICE-(row-1)/2-6]=pvertex;
      }
   }


   column_lim=LENGTH_FIRST_LATTICE;
   for (column=2 ; column<column_lim; column++){
      row=2*LENGTH_FIRST_LATTICE-1;
      vertex=vertex_from_coord(row, column);
      pvertex=web[vertex];
      pcell -> vertexlist[LENGTH_FIRST_LATTICE+column-4]=pvertex;

      row=1;
      vertex=vertex_from_coord(row, column);
      pvertex=web[vertex];
      pcell -> vertexlist[4*LENGTH_FIRST_LATTICE-column-7]=pvertex;
   }

   pcentercell=web_dual[vertex_from_coord(LENGTH_FIRST_LATTICE, LENGTH_FIRST_LATTICE/2)];
   if (!pcentercell){
      pcentercell=web_dual[vertex_from_coord(LENGTH_FIRST_LATTICE, LENGTH_FIRST_LATTICE/2) - 1];
      if (!pcentercell){
         locerror(" coord_gen ", " The central cell does not exist");
      }
   }
}
