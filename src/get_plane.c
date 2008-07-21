// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/************** get_plane() ******************/
/* this finds the plane defined by 3 coord_3D structs */
plane get_plane(coord_3D one, coord_3D two, coord_3D three){
plane p;
p.A=(two.j-one.j)*(three.k-one.k)-(three.j-one.j)*(two.k-one.k); 
p.B=(two.k-one.k)*(three.i-one.i)-(three.k-one.k)*(two.i-one.i); 
p.C=(two.i-one.i)*(three.j-one.j)-(three.i-one.i)*(two.j-one.j); 
p.D=-(p.A*one.i) - (p.B*one.j) - (p.C*one.k);
//dprint_coord_3D(one);
//dprint_coord_3D(two);
//dprint_coord_3D(three);
//dprint_plane(&p);
return p;
}
