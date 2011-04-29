// Function written by Oliver Grant, 2011
#include <glylib.h>
//#include <mylib.h>
//#include <molecules.h>
/************** get_distance_from_point_to_plane() ******************/
double get_distance_from_point_to_plane(plane p, coord_3D pt, int absl){

double sum, sqroot, d;

sum=(p.A)*(pt.i) + (p.B)*(pt.j) + (p.C)*(pt.k);
sqroot=sqrt(p.A*p.A + p.B*p.B + p.C*p.C);
d=sum/sqroot;

//get absolute value
if (absl==1){
	if (d<0) {d=d*(-1);}
	}
return d;
}
