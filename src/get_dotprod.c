// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* This function returns the dot product of two vectors. */
/*********************** get_dotprod() *******************/
double get_dotprod(vectormag_3D a,vectormag_3D b){
double dp=0;
dp=(a.i*b.i+a.j*b.j+a.k*b.k); 
return dp; 
}


