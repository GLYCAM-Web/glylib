// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* this turns a coordinate set into a vector */
vectormag_3D coord_to_vec(coord_3D c){
vectormag_3D nv;

nv.i=c.i;
nv.j=c.j;
nv.k=c.k;
// don't trust entry in 'd' -- normalize now to be certain
nv.d=sqrt(nv.i*nv.i+nv.j*nv.j+nv.k*nv.k);

return nv;
}

/* this turns a vector into a coordinate set */
coord_3D vec_to_coord(vectormag_3D v){
coord_3D nc;

nc.i=v.i;
nc.j=v.j;
nc.k=v.k;

return nc;
}
