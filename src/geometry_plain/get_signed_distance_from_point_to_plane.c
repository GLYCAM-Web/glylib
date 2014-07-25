// Function written by Spandana Makeneni, 2011
#include <glylib.h>
//#include <mylib.h>
//#include <molecules.h>
/************** get_distance_from_point_to_plane() ******************/
double get_signed_distance_from_point_to_plane(plane p, coord_3D pt ){

double sum, sqroot, d;
//printf("%f %f %f\n",pt.i,pt.j,pt.k);
//printf("%f %f %f\n",p.A,p.B,p.C);
//printf("%f %f %f\n",p.A*pt.i,p.B*pt.j,p.C*pt.k);
sum=(p.A)*(pt.i) + (p.B)*(pt.j) + (p.C)*(pt.k)+p.D;
//printf(" sum is %f\n",sum);
sqroot=sqrt(p.A*p.A + p.B*p.B + p.C*p.C);
d=sum/sqroot;

//get absolute value
return d;
}
