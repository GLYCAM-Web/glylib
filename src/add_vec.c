// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* this adds one vector to another (and sets d, too) */
vectormag_3D add_vec(vectormag_3D va,vectormag_3D vb){
vectormag_3D nv;

nv.i=va.i+vb.i;
nv.j=va.j+vb.j;
nv.k=va.k+vb.k;

// don't trust the d in nv
nv.d=sqrt(nv.i*nv.i+nv.j*nv.j+nv.k*nv.k);

return nv;
}

/* this subtracts vector vb from vector va (and sets d, too) */
vectormag_3D subtract_vec(vectormag_3D va,vectormag_3D vb){
vectormag_3D nv;

nv.i=va.i-vb.i;
nv.j=va.j-vb.j;
nv.k=va.k-vb.k;

// don't trust the d in nv
nv.d=sqrt(nv.i*nv.i+nv.j*nv.j+nv.k*nv.k);

return nv;
}
