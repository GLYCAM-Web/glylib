// Function written by B. Lachele Foley, 2007
// updated by same 20100307
//#include <mylib.h>
//#include <molecules.h>
#include "../inc/mylib.h"
#include "../inc/molecules.h"
/****************** get_residue_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Similar functions exist for the other
 * structures.
 */
coord_3D get_residue_COM(residue *r,atype *ATYPE, int xs){
int ta=0;
double tmm=0,*mass;
coord_3D c;
// Set the reference masses to use
// an calculate the mw while you're at it
mass=(double*)calloc(r[0].na,sizeof(double));
// calculate the molecular weight 
for(ta=0;ta<r[0].na;ta++){
	if((ATYPE!=NULL)&&(ATYPE!=0x0)){mass[ta]=ATYPE[r[0].a[ta].t].m;}
	else{ if(r[0].a[ta].m==0){
			fprintf(stderr,"Warning: Atom found with zero mass in get_residue_COM\n");}
		mass[ta]=r[0].a[ta].m; }
	tmm+=mass[ta];
	}
//printf("the molecular weight is %f\n",tmm);
if(tmm==0){mywhine("Found residue with zero mw.\n(Have types/masses been assigned?");}
r[0].m=tmm;

c.i=0;
c.j=0;
c.k=0; 
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
// calculate the center of mass
if(xs==-1){
for(ta=0;ta<r[0].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=r[0].a[ta].x.i*mass[ta];
        c.j+=r[0].a[ta].x.j*mass[ta];
        c.k+=r[0].a[ta].x.k*mass[ta];
        }
	}
else{
for(ta=0;ta<r[0].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=r[0].a[ta].xa[xs].i*mass[ta];
        c.j+=r[0].a[ta].xa[xs].j*mass[ta];
        c.k+=r[0].a[ta].xa[xs].k*mass[ta];
        }
	}
free(mass);
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
double tmm=0,**mass;
coord_3D c;
// Set the reference masses to use
// an calculate the mw while you're at it
mass=(double**)calloc(m[0].nr,sizeof(double*));
// calculate the molecular weight 
for(tr=0;tr<m[0].nr;tr++){
	mass[tr]=(double*)calloc(m[0].r[tr].na,sizeof(double));
for(ta=0;ta<m[0].r[tr].na;ta++){
	if((ATYPE!=NULL)&&(ATYPE!=0x0)){mass[tr][ta]=ATYPE[m[0].r[tr].a[ta].t].m;}
	else{ if(m[0].r[tr].a[ta].m==0){
			fprintf(stderr,"Warning: Atom found with zero mass in get_molecule_COM\n");}
		mass[tr][ta]=m[0].r[tr].a[ta].m; }
	tmm+=mass[tr][ta];
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
        c.i+=m[0].r[tr].a[ta].x.i*mass[tr][ta];
        c.j+=m[0].r[tr].a[ta].x.j*mass[tr][ta];
        c.k+=m[0].r[tr].a[ta].x.k*mass[tr][ta];
        }
	}
	}
else {
for(tr=0;tr<m[0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<m[0].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=m[0].r[tr].a[ta].xa[xs].i*mass[tr][ta];
        c.j+=m[0].r[tr].a[ta].xa[xs].j*mass[tr][ta];
        c.k+=m[0].r[tr].a[ta].xa[xs].k*mass[tr][ta];
        }
	}
	}
for(tr=0;tr<m[0].nr;tr++){
	free(mass[tr]);
	}
free(mass);
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
c.i/=tmm;
c.j/=tmm;
c.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return c;
}

/****************** get_assembly_molecule_COM() *********************/
/* Finds the center of mass based on the molecules and places it in 
 * the appropriate location.  Similar functions exist for the other
 * structures.
 */
coord_3D get_assembly_molecule_COM(assembly *a,atype *ATYPE, int xs){
int ta=0,tr=0,tm=0;
double tmm=0,***mass; 
coord_3D c;
// Set the reference masses to use
// and calculate the mw while you're at it
mass=(double***)calloc(a[0].nm,sizeof(double**));
for(tm=0;tm<a[0].nm;tm++){
	mass[tm]=(double**)calloc(a[0].m[tm][0].nr,sizeof(double*));
for(tr=0;tr<a[0].m[tm][0].nr;tr++){
	mass[tm][tr]=(double*)calloc(a[0].m[tm][0].r[tr].na,sizeof(double));
for(ta=0;ta<a[0].m[tm][0].r[tr].na;ta++){
// calculate the molecular weight 
	if((ATYPE!=NULL)&&(ATYPE!=0x0)){mass[tm][tr][ta]=ATYPE[a[0].m[tm][0].r[tr].a[ta].t].m;}
	else{ if(a[0].m[tm][0].r[tr].a[ta].m==0){
			fprintf(stderr,"Warning: Atom found with zero mass in get_assembly_molecule_COM\n");}
		mass[tm][tr][ta]=a[0].m[tm][0].r[tr].a[ta].m; }
	tmm+=mass[tm][tr][ta];
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
        c.i+=a[0].m[tm][0].r[tr].a[ta].x.i*mass[tm][tr][ta];
        c.j+=a[0].m[tm][0].r[tr].a[ta].x.j*mass[tm][tr][ta];
        c.k+=a[0].m[tm][0].r[tr].a[ta].x.k*mass[tm][tr][ta];
        }
	}
	}
	}
else {
for(tm=0;tm<a[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
for(tr=0;tr<a[0].m[tm][0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<a[0].m[tm][0].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=a[0].m[tm][0].r[tr].a[ta].xa[xs].i*mass[tm][tr][ta];
        c.j+=a[0].m[tm][0].r[tr].a[ta].xa[xs].j*mass[tm][tr][ta];
        c.k+=a[0].m[tm][0].r[tr].a[ta].xa[xs].k*mass[tm][tr][ta];
        }
	}
	}
	}
for(tm=0;tm<a[0].nm;tm++){
for(tr=0;tr<a[0].m[tm][0].nr;tr++){
	free(mass[tm][tr]);
	}
	free(mass[tm]);
	}
free(mass);
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
c.i/=tmm;
c.j/=tmm;
c.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return c;
}

/****************** get_ensemble_COM() *********************/
/* Finds the center of mass of the ensemble and places it in the
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
double tmm=0,***mass; 
coord_3D c;
// Set the reference masses to use
// and calculate the mw while you're at it
mass=(double***)calloc(e[0].nm,sizeof(double**));
for(tm=0;tm<e[0].nm;tm++){
	mass[tm]=(double**)calloc(e[0].m[tm].nr,sizeof(double*));
for(tr=0;tr<e[0].m[tm].nr;tr++){
	mass[tm][tr]=(double*)calloc(e[0].m[tm].r[tr].na,sizeof(double));
for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){
// calculate the molecular weight 
	if((ATYPE!=NULL)&&(ATYPE!=0x0)){mass[tm][tr][ta]=ATYPE[e[0].m[tm].r[tr].a[ta].t].m;}
	else{ if(e[0].m[tm].r[tr].a[ta].m==0){
			fprintf(stderr,"Warning: Atom found with zero mass in get_ensemble_COM\n");}
		mass[tm][tr][ta]=e[0].m[tm].r[tr].a[ta].m; }
	tmm+=mass[tm][tr][ta];
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
        c.i+=e[0].m[tm].r[tr].a[ta].x.i*mass[tm][tr][ta];
        c.j+=e[0].m[tm].r[tr].a[ta].x.j*mass[tm][tr][ta];
        c.k+=e[0].m[tm].r[tr].a[ta].x.k*mass[tm][tr][ta];
        }
	}
	}
	}
else {
for(tm=0;tm<e[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
for(tr=0;tr<e[0].m[tm].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){ // COM = sum(m_i * r_i) / sum(m_i)
        c.i+=e[0].m[tm].r[tr].a[ta].xa[xs].i*mass[tm][tr][ta];
        c.j+=e[0].m[tm].r[tr].a[ta].xa[xs].j*mass[tm][tr][ta];
        c.k+=e[0].m[tm].r[tr].a[ta].xa[xs].k*mass[tm][tr][ta];
        }
	}
	}
	}
for(tm=0;tm<e[0].nm;tm++){
for(tr=0;tr<e[0].m[tm].nr;tr++){
	free(mass[tm][tr]);
	}
	free(mass[tm]);
	}
free(mass);
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
c.i/=tmm;
c.j/=tmm;
c.k/=tmm; 
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return c;
} 
