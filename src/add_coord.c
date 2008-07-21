// Function written by B. Lachele Foley, 2007
#include <mylib.h>
//#include "../inc/mylib.h"
#include <molecules.h>
//#include "../inc/molecules.h"
/* this adds one coordinate to another */
coord_3D add_coord(coord_3D ca,coord_3D cb){
coord_3D nc;

nc.i=ca.i+cb.i;
nc.j=ca.j+cb.j;
nc.k=ca.k+cb.k;

return nc;
}

/* this subtracts vector vb from vector va (and sets d, too) */
coord_3D subtract_coord(coord_3D ca,coord_3D cb){
coord_3D nc;

nc.i=ca.i-cb.i;
nc.j=ca.j-cb.j;
nc.k=ca.k-cb.k;

return nc;
}

/* this one multiplies a coordinate by a scalar */
coord_3D scalarmult_coord(coord_3D ca,double cb){
coord_3D nc;

nc.i=ca.i*cb;
nc.j=ca.j*cb;
nc.k=ca.k*cb;

return nc;
}
