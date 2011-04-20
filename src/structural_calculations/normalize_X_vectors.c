// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* This function normalizes all the vectors in a molecule's atom
structures, in the location indicated by n1, and places the 
normalized forms in n2. The vectors are all normalized to one,
e.g., they are made into unit vectors. */
/************** normalize_molecule_vectors **********************/
void normalize_molecule_vectors(molecule *m,int n1, int n2){
int na=0,nb=0,nalloc=0;
int nuse=0,nsz=0,nbins=0;
// figure how big the vector array needs to be
if(n1>n2){nalloc=n1+1;}
else{nalloc=n2+1;} 
// normalize the vectors
for(na=0;na<m[0].nr;na++){ // for each residue
	for(nb=0;nb<m[0].r[na].na;nb++){ // for each atom
		//make sure there is space...
		nuse=malloc_usable_size(m[0].r[na].a[nb].v);
		nsz=sizeof(vectormag_3D);
		if(nsz==0){mywhine("memory issue (sizeof(vectormag_3D) given as zero) in normalize_molecule_vectors.c");}
		nbins=nuse/nsz;
		if(nbins<nalloc){mywhine("insufficient memory in normalize_molecule_vectors");}
		m[0].r[na].a[nb].v[n2]=normalize_vec(m[0].r[na].a[nb].v[n1]);
		}  // for each atom
	} // for each residue 
return;
} 
/************** normalize_ensemble_vectors **********************/
void normalize_ensemble_vectors(ensemble *e,int n1, int n2){
int nm=0;
// call for each molecule
for(nm=0;nm<e[0].nm;nm++){ normalize_molecule_vectors(&e[0].m[nm],n1,n2); } 
return;
} 
