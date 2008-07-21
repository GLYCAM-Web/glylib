// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
/* This function returns the cross product of two vectors.
  It also calculates and saves d in the structure. */
/******************** get_crossprod() *********************/
vectormag_3D get_crossprod(vectormag_3D a,vectormag_3D b){
vectormag_3D xp; 
xp.i=a.j*b.k-a.k*b.j;
xp.j=-(a.i*b.k-a.k*b.i);
xp.k=a.i*b.j-a.j*b.i;
xp.d=sqrt(xp.i*xp.i+xp.j*xp.j+xp.k*xp.k); 
return xp;
}
