// Function written by B. Lachele Foley, 2007
// ...and modified by BLFoley on 20071120 so that 
//   the functions do not require arguments
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* this zero's a vector */
vectormag_3D zero_vec(){
vectormag_3D nv;

nv.i=nv.j=nv.k=nv.d=0;

return nv;
}


/* this zero's a coordinate set */
coord_3D zero_coord(){
coord_3D nc;

nc.i=nc.j=nc.k=0;

return nc;
}
