// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/************** get_ring_center() ******************/
/* this finds the center point for a set of atoms.
 c is a pointer to a list of nc coord_3D sets */
coord_3D  get_geometric_center(coord_3D *c,int nc){
int gca=0;
coord_3D C;

if(nc==0){mywhine("cannot get center for zero points (get_geometric_center)");}
C.i=C.j=C.k=0;
//printf("initialized:\n");
//dprint_coord_3D(&C);
//printf("in loop:\n");
for(gca=0;gca<nc;gca++){
//      X.i (center) = avg (atomscoords.x.i) ; and so forth
	C.i+=c[gca].i; 
	C.j+=c[gca].j; 
	C.k+=c[gca].k; 
	//dprint_coord_3D(&C);
	}
C.i/=nc;
C.j/=nc;
C.k/=nc;
//printf("done:\n");
//dprint_coord_3D(&C);
return C;
}
/* This is just like get_geometric_center except it expects an array
 of double pointers.  Uses less memory when referring to coordinates
 within the larger structures.  Each pointer in the array should 
 point to a coordinate to be used in the calculation. */
coord_3D  get_geometric_center_dp(coord_3D **c,int nc){
int gca=0;
coord_3D C;

if(nc==0){mywhine("cannot get center for zero points (get_geometric_center)");}
C.i=C.j=C.k=0;
//printf("initialized:\n");
//dprint_coord_3D(&C);
//printf("in loop:\n");
for(gca=0;gca<nc;gca++){
//      X.i (center) = avg (atomscoords.x.i) ; and so forth
	C.i+=c[gca][0].i; 
	C.j+=c[gca][0].j; 
	C.k+=c[gca][0].k; 
	//dprint_coord_3D(&C);
	}
C.i/=nc;
C.j/=nc;
C.k/=nc;
//printf("done:\n");
//dprint_coord_3D(&C);
return C;
}

