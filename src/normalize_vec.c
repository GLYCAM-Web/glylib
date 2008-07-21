// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* this normalizes a vector */
vectormag_3D normalize_vec(vectormag_3D v){
vectormag_3D nv;
double temp=0;

// don't trust entry in 'd' -- normalize now to be certain
temp=sqrt(v.i*v.i+v.j*v.j+v.k*v.k);

if(temp!=0){
	nv.i=v.i/temp;
	nv.j=v.j/temp;
	nv.k=v.k/temp;
	nv.d=1;
	}
else{ nv.i=nv.j=nv.k=nv.d=0; }

return nv;
}
