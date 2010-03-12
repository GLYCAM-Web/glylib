// Function written by B. Lachele Foley, 2007
//
// Deprecated.  Use "set_COM" functions instead
//
//#include <mylib.h>
//#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/****************** assign_residue_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Similar functions exist for the other
 * structures.
 */
/*
void assign_residue_COM(residue *r,atype *ATYPE){
int ta=0;
double tmm=0; 
// calculate the molecular weight 
for(ta=0;ta<r[0].na;ta++){
	tmm+=ATYPE[r[0].a[ta].t].m;
	}
//printf("the molecular weight is %f\n",tmm);
if(tmm==0){mywhine("molecular weight of molecule %s is zero.\n(Have types/masses been assigned?");}
r[0].m=tmm;

r[0].COM.i=0;
r[0].COM.j=0;
r[0].COM.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
for(ta=0;ta<r[0].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        r[0].COM.i+=r[0].a[ta].x.i*ATYPE[r[0].a[ta].t].m;
        r[0].COM.j+=r[0].a[ta].x.j*ATYPE[r[0].a[ta].t].m;
        r[0].COM.k+=r[0].a[ta].x.k*ATYPE[r[0].a[ta].t].m;
        }
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
r[0].COM.i/=tmm;
r[0].COM.j/=tmm;
r[0].COM.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return;
}
*/

/****************** assign_molecule_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Similar functions exist for the other
 * structures.
 */
/* void assign_molecule_COM(molecule *m,atype *ATYPE){
int ta=0,tr=0;
double tmm=0; 
// calculate the molecular weight 
for(tr=0;tr<m[0].nr;tr++){
for(ta=0;ta<m[0].r[tr].na;ta++){
	tmm+=ATYPE[m[0].r[tr].a[ta].t].m;
	}
	}
//printf("the molecular weight is %f\n",tmm);
if(tmm==0){mywhine("molecular weight of molecule %s is zero.\n(Have types/masses been assigned?");}
m[0].m=tmm;

m[0].COM.i=0;
m[0].COM.j=0;
m[0].COM.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
for(tr=0;tr<m[0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<m[0].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        m[0].COM.i+=m[0].r[tr].a[ta].x.i*ATYPE[m[0].r[tr].a[ta].t].m;
        m[0].COM.j+=m[0].r[tr].a[ta].x.j*ATYPE[m[0].r[tr].a[ta].t].m;
        m[0].COM.k+=m[0].r[tr].a[ta].x.k*ATYPE[m[0].r[tr].a[ta].t].m;
        }
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
m[0].COM.i/=tmm;
m[0].COM.j/=tmm;
m[0].COM.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return;
}
*/

/****************** assign_assembly_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Similar functions exist for the other
 * structures.
 */
/*void assign_assembly_COM(assembly *a,atype *ATYPE){
int ta=0,tr=0,tm=0;
double tmm=0; 
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

a[0].COM.i=0;
a[0].COM.j=0;
a[0].COM.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
for(tm=0;tm<a[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
for(tr=0;tr<a[0].m[tm][0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<a[0].m[tm][0].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        a[0].COM.i+=a[0].m[tm][0].r[tr].a[ta].x.i*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        a[0].COM.j+=a[0].m[tm][0].r[tr].a[ta].x.j*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        a[0].COM.k+=a[0].m[tm][0].r[tr].a[ta].x.k*ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;
        }
	}
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
a[0].COM.i/=tmm;
a[0].COM.j/=tmm;
a[0].COM.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return;
}
*/
/****************** assign_ensemble_COM() *********************/
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
/* void assign_ensemble_COM(ensemble *e,atype *ATYPE){
int ta=0,tr=0,tm=0;
double tmm=0; 
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

e[0].COM.i=0;
e[0].COM.j=0;
e[0].COM.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
for(tm=0;tm<e[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
for(tr=0;tr<e[0].m[tm].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        e[0].COM.i+=e[0].m[tm].r[tr].a[ta].x.i*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        e[0].COM.j+=e[0].m[tm].r[tr].a[ta].x.j*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        e[0].COM.k+=e[0].m[tm].r[tr].a[ta].x.k*ATYPE[e[0].m[tm].r[tr].a[ta].t].m;
        }
	}
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
e[0].COM.i/=tmm;
e[0].COM.j/=tmm;
e[0].COM.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return;
} 
*/
