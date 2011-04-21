// Function written by Oliver Grant, 2011
#include <mylib.h>
#include <molecules.h>
/************** get_angle_between_plane_and_vector() ******************/
double get_angle_between_plane_and_vector(plane p, coord_3D pt1, coord_3D pt2){

double magnitude, angle;
vectormag_3D plUV,ptV,ptUV; //plane unit vector, points vector, points unit vector

//get normal unit vector of the plane
magnitude=sqrt(p.A*p.A + p.B*p.B + p.C*p.C); //can i trust p.D?
plUV.i=p.A/magnitude;
plUV.j=p.B/magnitude;
plUV.k=p.C/magnitude;

//get vector of the two points
ptV.i=(pt1.i-pt2.i);
ptV.j=(pt1.j-pt2.j);
ptV.k=(pt1.k-pt2.k);
magnitude=sqrt((ptV.i*ptV.i) + (ptV.j*ptV.j) + (ptV.k*ptV.k));

// normalize this vector
ptUV.i=ptV.i/magnitude;
ptUV.j=ptV.j/magnitude;
ptUV.k=ptV.k/magnitude;

printf("magnitude=%f\nptUV.i=%f,ptUV.j=%f,ptUV.k=%f,plUV.i=%f,plUV.j=%f,plUV.k=%f\n",magnitude,ptUV.i,ptUV.j,ptUV.k,plUV.i,plUV.j,plUV.k);

// orientation of points vector with respect to the plane
angle=(get_dotprod(plUV,ptUV));  /* angle in radians */

return angle;
}
