// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/****************** translate_to_COM() *********************/
/* Finds the COM for molecule m, saves the COM coordinates
in the structure and translates the molecule so that the xyz
origin is at the COM.  The transformed coordinates can go in
the original coordinate spot (tw=-1) or in one of the alternate
coordinate spots (tw>=0). */

// START HERE -- retire this function ASAP (rewrite others not to use)

void translate_to_COM(molecule *m,atype *ATYPE,int tw){
int ta=0,tuse=0,tsz=0,tbins=0;
double tmm=0; 

// calculate the molecular weight of residue 0 (the molecule, here) 
for(ta=0;ta<m[0].r[0].na;ta++){
	tmm+=ATYPE[m[0].r[0].a[ta].t].m;
	}
//printf("the molecular weight is %f\n",tmm);
if(tmm==0){mywhine("molecular weight of molecule %s is zero.\n(Have types/masses been assigned?");}
// glucose $molwt=174.078;

m[0].COM.i=0;
m[0].COM.j=0;
m[0].COM.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
for(ta=0;ta<m[0].r[0].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        m[0].COM.i+=m[0].r[0].a[ta].x.i*ATYPE[m[0].r[0].a[ta].t].m;
        m[0].COM.j+=m[0].r[0].a[ta].x.j*ATYPE[m[0].r[0].a[ta].t].m;
        m[0].COM.k+=m[0].r[0].a[ta].x.k*ATYPE[m[0].r[0].a[ta].t].m;
        }
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
m[0].COM.i/=tmm;
m[0].COM.j/=tmm;
m[0].COM.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);

// move center of mass to specified location
if(tw==-1){
	for(ta=0;ta<m[0].r[0].na;ta++){ // xyz(new) = xyz(old) - COM
        	m[0].r[0].a[ta].x.i-=m[0].COM.i;
        	m[0].r[0].a[ta].x.j-=m[0].COM.j;
        	m[0].r[0].a[ta].x.k-=m[0].COM.k;
        	}
	}
else{
	for(ta=0;ta<m[0].r[0].na;ta++){ // xyz(new) = xyz(old) - COM
		tuse=malloc_usable_size(m[0].r[0].a[ta].xa);
                tsz=sizeof(coord_3D);
		if(tsz==0){mywhine("memory issue: sizeof(coord_3D) given as zero in translate_to_COM.c");}
                tbins=tuse/tsz; 
		//print_coord_3D(&m[0].r[0].a[ta].x);
		if(tbins<(tw+1)){mywhine("insufficient memory allocated before translate_to_COM");}
        	m[0].r[0].a[ta].xa[tw].i=m[0].r[0].a[ta].x.i-m[0].COM.i;
        	m[0].r[0].a[ta].xa[tw].j=m[0].r[0].a[ta].x.j-m[0].COM.j;
        	m[0].r[0].a[ta].xa[tw].k=m[0].r[0].a[ta].x.k-m[0].COM.k;
        	}
	}
return;
}
