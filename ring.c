#include <mylib.h>
#include <molecules.h>
#include <math.h>
plane get_plane_for_ring(int n,coord_3D *r);

int main(int argc,char *argv[]){
int n=6;
coord_3D *m;
m=(coord_3D*)calloc(n,sizeof(coord_3D));//array of coordinates 
m[0].i = 0.1607776;
m[0].j = 1.4429037;
m[0].k = 1.5216059;
m[1].i = 0.1596404;
m[1].j = 1.4967956;
m[1].k = 3.04689463;
m[2].i = -1.219758;
m[2].j = 1.077754248;
m[2].k = 3.545446357;
m[3].i = -1.466034;
m[3].j = -0.342773;
m[3].k = 3.0468395;
m[4].i = -1.412516;
m[4].j = -0.337318;
m[4].k = 1.5216409;
m[5].i = -0.128463;
m[5].j = 0.1132565;
m[5].k = 1.0682333;
plane newp;
newp = get_plane_for_ring(n,m);
printf("%f %f %f %f\n",newp.A,newp.B,newp.C,newp.D);
return 0;

}
