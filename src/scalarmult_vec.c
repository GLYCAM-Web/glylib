// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
/* this multiplies a vector by a scalar (and sets d, too) */
vectormag_3D scalarmult_vec(vectormag_3D v,double s){
vectormag_3D nv;

nv.i=v.i*s;
nv.j=v.j*s;
nv.k=v.k*s;

// don't trust the d in nv
nv.d=sqrt(nv.i*nv.i+nv.j*nv.j+nv.k*nv.k);

return nv;
}
