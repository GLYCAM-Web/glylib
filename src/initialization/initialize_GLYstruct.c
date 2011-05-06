/** \file initialize_GLYstruct.c 
 * Provides initialization functions for structures in molecules.h
 *
 * Probably needs some serious updating (20080813, BLF).
 *
 * NOTES:	char* variables are not initialized
 *		char** are treated like "all other pointers" (below)
 *		void pointers are not initialized
 *		all other pointers:
 *			-- the number of pointers integer is set to zero
 *			-- one pointer is calloc'd (in anticipation of realloc)
 * 
 * begun on 20071207 by BLFoley 
 */
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"

#define INDEX_INIT -100000000 ///< To catch use of un-set indices early on.

void initialize_coord_3D(coord_3D *c) {
c[0].i=c[0].j=c[0].k=0; 
return;
}

void initialize_vectormag_3D(vectormag_3D *v) {
v[0].i=v[0].j=v[0].k=v[0].d=0; 
return;
}

void initialize_plane(plane *p) {
p[0].A=p[0].B=p[0].C=p[0].D=0;
return;
}

void initialize_atype(atype *at) {
at[0].n=0; // atomic number
at[0].nb=0; // number of typical bonds for this atom type
at[0].nlp=0; // number of typical lone pairs (to check geometry sanity)
at[0].m=0; // mass of element type
at[0].bo=(double*)calloc(1,sizeof(double)); // typical bond orders (nb of these)
return;
} 

void initialize_rtype(rtype *rt) {
rt[0].c=0; // main class (amino acid, glycan, solvent, etc.)
rt[0].nac=0;// number alternate classes
rt[0].ac=(int*)calloc(1,sizeof(int)); // number alternate classes, those classes
rt[0].nVP=0; // number of other information
return;
}

void initialize_mtype(mtype *mt) {
mt[0].c=0; // type class (amino acid, glycan, solvent, etc.)
mt[0].nVP=0; // number of other information
return;
}

void initialize_types(types *t) {
t[0].na=0; // # of atom types
t[0].a=(atype*)calloc(1,sizeof(atype)); // na of these
t[0].nr=0; // # of residue types
t[0].r=(rtype*)calloc(1,sizeof(rtype)); // nr of these
t[0].nm=0; // # of molecule types
t[0].m=(mtype*)calloc(1,sizeof(mtype)); // nm of these
return;
}

void initialize_molindex(molindex *mi) {
mi[0].i=INDEX_INIT; // general index
mi[0].m=INDEX_INIT; // molecule index
mi[0].r=INDEX_INIT; // residue index
mi[0].a=INDEX_INIT; // atom index
return;
}

void initialize_ensindex(ensindex *ei) {
ei[0].i=INDEX_INIT; // general index
ei[0].E=INDEX_INIT; // ensemble
ei[0].A=INDEX_INIT; // assembly
ei[0].m=INDEX_INIT; // molecule index
ei[0].r=INDEX_INIT; // residue index
ei[0].a=INDEX_INIT; // atom index
return;
}

void initialize_bond(bond *b) {
b[0].s.m=INDEX_INIT; // "source" -- index to first atom in bond
b[0].s.r=INDEX_INIT; // "source" -- index to first atom in bond
b[0].s.a=INDEX_INIT; // "source" -- index to first atom in bond
b[0].t.m=INDEX_INIT; // "target" -- index to the other atom in the bond
b[0].t.r=INDEX_INIT; // "target" -- index to the other atom in the bond
b[0].t.a=INDEX_INIT; // "target" -- index to the other atom in the bond
b[0].o=0; // order of bond
return;
}

void initialize_bondset(bondset *bs) {
bs[0].n=0; // number of bonds
bs[0].b=(bond*)calloc(1,sizeof(bond)); // n of these
return;
}

void initialize_molbond(molbond *mb) {
mb[0].s.i=mb[0].s.m=mb[0].s.r=mb[0].s.a=INDEX_INIT;
mb[0].t.i=mb[0].t.m=mb[0].t.r=mb[0].t.a=INDEX_INIT;
mb[0].o=0; // order
return;
}

void initialize_molbondset(molbondset *mbs) {
mbs[0].n=0;  // source & target atoms for bond
mbs[0].b=(molbond*)calloc(1,sizeof(molbond)); // order
return;
}

void initialize_atom(atom *a) {
a[0].n=0; // atom number or other identifying index
a[0].N=NULL;
a[0].T=NULL;
a[0].D=NULL;
a[0].E=NULL; 
a[0].cID=NULL; 
a[0].typ=NULL; 
a[0].m=0;
a[0].t=0; // type number -- must correspond to assignments of "atype" (see) 
initialize_molindex(&a[0].moli);
a[0].nb=0; // number of actual bonds (not expected bonds)
a[0].b=(bond*)calloc(1,sizeof(bond)); // bond structures (nb of these)
a[0].nmb=0; // number of bonds to other residues or molecules
a[0].mb=(molbond*)calloc(1,sizeof(molbond)); // nmb of these
a[0].mTi=-1;
a[0].rTi=-1;
a[0].x.i=a[0].x.j=a[0].x.k=0; // atom's coordinates 
a[0].xv.i=a[0].xv.j=a[0].xv.k=0; // atom's coordinates 
a[0].nalt=0; // number of alternate coordinate sets
a[0].xa=(coord_3D*)calloc(1,sizeof(coord_3D)); // nalt of alternate coords
a[0].nxva=0; // number of alternate coordinate sets
a[0].xva=(coord_3D*)calloc(1,sizeof(coord_3D)); // nalt of alternate coords
a[0].nvec=0; // number of vector sets
a[0].v=(vectormag_3D*)calloc(1,sizeof(vectormag_3D)); // vector sets
a[0].nch=0;
a[0].ch=(double*)calloc(1,sizeof(double));
a[0].ni=0; // number of other indices
a[0].i=(int*)calloc(1,sizeof(int)); // other indices, as needed (ni of these)
a[0].nd=0; // number of double-precision parameters
a[0].d=(double*)calloc(1,sizeof(double)); // other parameters, as needed (nd of these) 
a[0].nensi=0; // number of double-precision parameters
a[0].ensi=(ensindex*)calloc(1,sizeof(ensindex)); // other parameters, as needed (nd of these)
a[0].nOD=0; // number of double-precision parameters
a[0].OD=(char**)calloc(1,sizeof(char*)); // other parameters, as needed (nd of these) 
a[0].nVP=0; // number of void structures
return;
}

// START HERE -- UPDATE all of the residue, molecule and higher structures.
void initialize_residue(residue *r) {
r[0].n=0; // residue number given in input file
r[0].cID=NULL;
r[0].IC=NULL;
r[0].N=NULL;
r[0].T=NULL;
r[0].D=NULL;
r[0].altname=NULL;
r[0].typ=NULL;
r[0].t=0; // index for rtype
initialize_molindex(&r[0].moli);
r[0].na=0; // number of atoms in residue
r[0].m=0; // molecular weight
r[0].COM.i=r[0].COM.j=r[0].COM.k=0; // center of mass for molecule
r[0].a=(atom*)calloc(1,sizeof(atom)); // atom structures (na of these)
r[0].aT=(atom_node*)calloc(1,sizeof(atom_node)); // atom structures (na of these)
r[0].nrb=0; // number of bonds to other residues or molecules
r[0].rb=(molbond*)calloc(1,sizeof(molbond)); // nmb of these
r[0].mTi=-1;
r[0].nbs=0; // number of bond sets 
r[0].bs=(molbondset*)calloc(1,sizeof(molbondset)); // (consecutive bonds, use these for plotting, etc.)
r[0].nring=0; // number of simple rings (no cage structures, etc.)
r[0].nrc=0; // number of ring/reference coordinate sets defined
r[0].rc=(coord_3D*)calloc(1,sizeof(coord_3D)); // coordinates for ring/reference centers
r[0].nrp=0; // number of ring planes defined
r[0].rp=(plane*)calloc(1,sizeof(plane)); // equations for average/approximate/exact/etc. ring planes (where useful)
r[0].ni=0; // number of other indices
r[0].i=(int*)calloc(1,sizeof(int)); // other indices, as needed (ni of these)
r[0].nd=0; // number of double-precision parameters
r[0].d=(double*)calloc(1,sizeof(double)); // other parameters, as needed (nd of these)
r[0].nOD=0; // number of double-precision parameters
r[0].OD=(char**)calloc(1,sizeof(char*)); // other parameters, as needed (nd of these)
r[0].nVP=0; // number of void structures
return;
}

void initialize_molecule(molecule *m) {
m[0].i=0; // index
m[0].t=0; // index
m[0].m=0; // molecular weight
m[0].COM.i=m[0].COM.j=m[0].COM.k=0; // center of mass for molecule
m[0].na=0; // total number of atoms in molecule
m[0].nr=0; // number of residues
m[0].r=(residue*)calloc(1,sizeof(residue)); // pointers to residues
/*
m[0].nrb=0; // number of bonds between residues
m[0].rb=(molbond*)calloc(1,sizeof(molbond)); // nrb of these descriptions of bonds 
*/
m[0].nrbs=0; // number of sets of bonds between residues (for example, linear chains)
m[0].rbs=(molbondset*)calloc(1,sizeof(molbondset)); // nrbs of these sets
m[0].nrc=0; // number of additional reference coordinates (rings, for example)
m[0].rc=(coord_3D*)calloc(1,sizeof(coord_3D)); // nrc of these
m[0].nBOX=0; // changed to new BOX member on 20080813 BLF
m[0].noi=0; // number of other indices
m[0].oi=(int*)calloc(1,sizeof(int)); // other indices, as needed (ni of these)
m[0].nd=0; // number of double-precision parameters
m[0].d=(double*)calloc(1,sizeof(double)); // other parameters, as needed (nd of these)
m[0].nVP=0; // number of void structures
return;
}

void initialize_dockinfo(dockinfo *di){
printf("This program will exit now, because it requests initialization of\n");
printf("a dockinfo structure, but dockinfo structures should be initialized\n");
printf("within the main program.\n");
exit(1);
return;
}

void initialize_assembly(assembly *A) { // structure for groups of molecules within a larger structure
A[0].i=0; // index
A[0].mass=0; // mass of assembly
A[0].COM.i=A[0].COM.j=A[0].COM.k=0; // center of mass 
A[0].nm=0; // number of molecule structures
A[0].m=(molecule**)calloc(1,sizeof(molecule*)); // nm of these
A[0].nb=0; // number of bonds/connections between molecules (H-bonds, for example)
A[0].b=(molbond*)calloc(1,sizeof(molbond)); // nmb of these descriptions of connection 
A[0].nmbs=0; // number of sets of connections between molecules (for example, linear chains)
A[0].mbs=(molbondset*)calloc(1,sizeof(molbondset)); // nmbs of these sets
A[0].nBOX=0; ///< changed to new BOX member on 20080813 BLF
//A[0].boxl.i=A[0].boxl.j=A[0].boxl.k=0; 
//A[0].boxh.i=A[0].boxh.j=A[0].boxh.k=0; 
A[0].nVP=0; // number of void structures
return;
}

void initialize_ensemble(ensemble *E) { // structure for a larger group of molecules, assemblies, etc.
E[0].i=0; // index
E[0].mass=0; // mass of ensemble
E[0].COM.i=E[0].COM.j=E[0].COM.k=0; // center of mass 
E[0].nm=0; // number of molecule structures
E[0].m=(molecule*)calloc(1,sizeof(molecule)); // nm of these
E[0].nA=0; // number of assembly structures
E[0].A=(assembly**)calloc(1,sizeof(assembly*)); // na of these
E[0].nBOX=0; ///< changed to new BOX member on 20080813 BLF
//E[0].boxl.i=E[0].boxl.j=E[0].boxl.k=0; // center of mass 
//E[0].boxh.i=E[0].boxh.j=E[0].boxh.k=0; // center of mass 
E[0].nVP=0; // number of void structures
return;
}

