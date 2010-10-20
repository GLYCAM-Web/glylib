/***Function written By Spandana Makeneni 2010***/
#include<mylib.h>
#include<molecules.h>
/***This defines the average plane for an array of coordinates according to puckering paramaters***/

plane get_plane_for_ring(int n,coord_3D *r){
int jval[n];
int l;
plane pval;//plane for returning the values of plane 
vectormag_3D Rj,Rcos,Rsin,R1,R2,R1xR2,avg_coord;
R1.i=R1.j=R1.k=R2.i=R2.j=R2.k=0;
avg_coord.i=avg_coord.j=avg_coord.k=0;
for(l=0;l<n;l++){
	jval[l]=l+1;//getting j vals for further equations
	}
for(l=0;l<n;l++){//for loop for calculating the Rj vals
	Rj.i= r[l].i;
	Rj.j= r[l].j;
	Rj.k= r[l].k;
	Rsin.i = (Rj.i*sin(2*PI*(jval[l]-1)/n));
	Rsin.j = (Rj.j*sin(2*PI*(jval[l]-1)/n));
	Rsin.k = (Rj.k*sin(2*PI*(jval[l]-1)/n));
	Rcos.i = (Rj.i*cos(2*PI*(jval[l]-1)/n));
	Rcos.j = (Rj.j*cos(2*PI*(jval[l]-1)/n));
	Rcos.k = (Rj.k*cos(2*PI*(jval[l]-1)/n));
	R1.i=R1.i+Rsin.i;
	R1.j=R1.j+Rsin.j;
	R1.k=R1.k+Rsin.k;
	R2.i=R2.i+Rcos.i;
	R2.j=R2.j+Rcos.j;
	R2.k=R2.k+Rcos.k;
	}
R1xR2 = get_crossprod(R1,R2);//cross prod of Rjs
printf("the R1xR2 is %f %f %f %f\n",R1xR2.i,R1xR2.j,R1xR2.k,R1xR2.d);
printf("the R1val is %f  %f\n",R1.i,R2.i);
printf("the R2val is %f  %f\n",R1.j,R2.j);
printf("the R3val is %f  %f\n",R1.k,R2.k);
pval.A = (R1xR2.i/R1xR2.d);
pval.B = (R1xR2.j/R1xR2.d);
pval.C = (R1xR2.k/R1xR2.d);
for(l=0;l<n;l++){//for loop for calculating the avg x,y,z coordinates
	avg_coord.i=avg_coord.i+r[l].i;
	avg_coord.j=avg_coord.j+r[l].j;
	avg_coord.k=avg_coord.k+r[l].k;
	}
avg_coord.i=(avg_coord.i/n);
avg_coord.j=(avg_coord.j/n);
avg_coord.k=(avg_coord.k/n);
printf("the avg coords are %f %f %f\n",avg_coord.i,avg_coord.j,avg_coord.k);
pval.D= -(pval.A*avg_coord.i+pval.B*avg_coord.j+pval.C*avg_coord.k);
printf("the dval is %f\n",pval.D);

return pval;

}
