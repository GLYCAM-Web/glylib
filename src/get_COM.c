// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/****************** get_residue_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Similar functions exist for the other
 * structures.
 */
coord_3D get_residue_COM(residue *r,atype *ATYPE, int xs){
int ta=0;
double tmm=0; 
coord_3D c;
// calculate the molecular weight 
for(ta=0;ta<r[0].na;ta++){
	tmm+=ATYPE[r[0].a[ta].t].m;
	}
//printf("the molecular weight is %f\n",tmm);
if(tmm==0){mywhine("molecular weight of molecule %s is zero.\n(Have types/masses been assigned?");}
r[0].m=tmm;

c.i=0;
c.j=0;
c.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
if(xs==-1){
for(ta=0;ta<r[0].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=r[0].a[ta].x.i*ATYPE[r[0].a[ta].t].m;
        c.j+=r[0].a[ta].x.j*ATYPE[r[0].a[ta].t].m;
        c.k+=r[0].a[ta].x.k*ATYPE[r[0].a[ta].t].m;
        }
	}
else{
for(ta=0;ta<r[0].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=r[0].a[ta].xa[xs].i*ATYPE[r[0].a[ta].t].m;
        c.j+=r[0].a[ta].xa[xs].j*ATYPE[r[0].a[ta].t].m;
        c.k+=r[0].a[ta].xa[xs].k*ATYPE[r[0].a[ta].t].m;
        }
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
c.i/=tmm;
c.j/=tmm;
c.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return c;
}

/****************** get_molecule_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Similar functions exist for the other
 * structures.
 */
coord_3D get_molecule_COM(molecule *m,atype *ATYPE, int xs){
int ta=0,tr=0;
double tmm=0; 
coord_3D c;
// calculate the molecular weight 
for(tr=0;tr<m[0].nr;tr++){
for(ta=0;ta<m[0].r[tr].na;ta++){
	tmm+=ATYPE[m[0].r[tr].a[ta].t].m;
	}
	}
//printf("the molecular weight is %f\n",tmm);
if(tmm==0){mywhine("molecular weight of molecule %s is zero.\n(Have types/masses been assigned?");}
m[0].m=tmm;

c.i=0;
c.j=0;
c.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
if(xs==-1){
for(tr=0;tr<m[0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<m[0].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=m[0].r[tr].a[ta].x.i*ATYPE[m[0].r[tr].a[ta].t].m;
        c.j+=m[0].r[tr].a[ta].x.j*ATYPE[m[0].r[tr].a[ta].t].m;
        c.k+=m[0].r[tr].a[ta].x.k*ATYPE[m[0].r[tr].a[ta].t].m;
        }
	}
	}
else {
for(tr=0;tr<m[0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<m[0].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=m[0].r[tr].a[ta].xa[xs].i*ATYPE[m[0].r[tr].a[ta].t].m;
        c.j+=m[0].r[tr].a[ta].xa[xs].j*ATYPE[m[0].r[tr].a[ta].t].m;
        c.k+=m[0].r[tr].a[ta].xa[xs].k*ATYPE[m[0].r[tr].a[ta].t].m;
        }
	}
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
c.i/=tmm;
c.j/=tmm;
c.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return c;
}

/****************** get_assembly_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Similar functions exist for the other
 * structures.
 */
coord_3D get_assembly_COM(assembly *a,atype *ATYPE, int xs){
int ta=0,tr=0,tm=0;
double tmm=0; 
coord_3D c;
// calculate the molecular weight 
for(tm=0;tm<a[0].nm;tm++){
for(tr=0;tr<a[0].m[tm][0].nr;tr++){
for(ta=0;ta<a[0].m[tm][0].r[tr].na;ta++){
	tmm+=ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
	}
	}
	}
//printf("the molecular weight is %f\n",tmm);
if(tmm==0){mywhine("molecular weight of molecule %s is zero.\n(Have types/masses been assigned?");}
a[0].mass=tmm;

c.i=0;
c.j=0;
c.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
if(xs==-1){
for(tm=0;tm<a[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
for(tr=0;tr<a[0].m[tm][0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<a[0].m[tm][0].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=a[0].m[tm][0].r[tr].a[ta].x.i*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        c.j+=a[0].m[tm][0].r[tr].a[ta].x.j*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        c.k+=a[0].m[tm][0].r[tr].a[ta].x.k*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        }
	}
	}
	}
else {
for(tm=0;tm<a[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
for(tr=0;tr<a[0].m[tm][0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<a[0].m[tm][0].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=a[0].m[tm][0].r[tr].a[ta].xa[xs].i*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        c.j+=a[0].m[tm][0].r[tr].a[ta].xa[xs].j*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        c.k+=a[0].m[tm][0].r[tr].a[ta].xa[xs].k*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        }
	}
	}
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
c.i/=tmm;
c.j/=tmm;
c.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return c;
}

/****************** get_ensemble_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Similar functions exist for the other
 * structures.
 * NOTE!!  Assemblies within ensembles are assumed to be redundant.
 * 	That is, the molecules in the ensemble are also represented
 * 	within the assemblies (copying the molecule pointers keeps 
 * 	the data from taking extra space).  So, when calculating the
 * 	center of mass for the ensemble, this function ONLY considers
 * 	molecules and not assemblies.
 */
coord_3D get_ensemble_COM(ensemble *e,atype *ATYPE, int xs){
int ta=0,tr=0,tm=0;
double tmm=0; 
coord_3D c;
// calculate the molecular weight 
for(tm=0;tm<e[0].nm;tm++){
for(tr=0;tr<e[0].m[tm].nr;tr++){
for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){
	tmm+=ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
	}
	}
	}
//printf("the molecular weight is %f\n",tmm);
if(tmm==0){mywhine("molecular weight of molecule %s is zero.\n(Have types/masses been assigned?");}
e[0].mass=tmm;

c.i=0;
c.j=0;
c.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
if(xs==-1){
for(tm=0;tm<e[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
for(tr=0;tr<e[0].m[tm].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=e[0].m[tm].r[tr].a[ta].x.i*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        c.j+=e[0].m[tm].r[tr].a[ta].x.j*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        c.k+=e[0].m[tm].r[tr].a[ta].x.k*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        }
	}
	}
	}
else {
for(tm=0;tm<e[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
for(tr=0;tr<e[0].m[tm].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=e[0].m[tm].r[tr].a[ta].xa[xs].i*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        c.j+=e[0].m[tm].r[tr].a[ta].xa[xs].j*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        c.k+=e[0].m[tm].r[tr].a[ta].xa[xs].k*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        }
	}
	}
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
c.i/=tmm;
c.j/=tmm;
c.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return c;
} 
