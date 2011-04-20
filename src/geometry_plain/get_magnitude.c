// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* this returns a vector with magnitude in d*/
double get_magnitude(vectormag_3D v){
double d=0;
d=sqrt(v.i*v.i+v.j*v.j+v.k*v.k); 
return d;
}
