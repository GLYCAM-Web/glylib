// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/****************** translate_by_XYZ() *********************/
/* Adds the vector vec to all coordinates in the specified
 * structure.  Uses xs (source) as the coordinates to move and
 * xd (destination) as the location for the new coordinates.
 * As usual xs or xd = -1 means main set.
 */
//
//The following are useful for troubleshooting memory issues
//tuse=malloc_usable_size(m[0].r[0].a[ta].xa);
//tsz=sizeof(coord_3D);
//if(tsz==0){mywhine("memory issue: sizeof(coord_3D) given as zero in translate_to_COM.c");}
//tbins=tuse/tsz; 
//print_coord_3D(&m[0].r[0].a[ta].x);
//if(tbins<(tw+1)){mywhine("insufficient memory allocated before translate_to_COM");}

void translate_residue_by_XYZ(residue *r,int xs,int xd,coord_3D vec){
int ta=0;
if(r[0].na<=0){mywhine("Residue contains no atoms in translate_residue_by_XYZ.");}
if((xs<-1)||(xd<-1)){mywhine("Invalid source or destination coordinates to translate_residue_by_XYZ.");}
if((xs==-1)&&(xd==-1)){
	for(ta=0;ta<r[0].na;ta++){ 
        	r[0].a[ta].x.i+=vec.i;
        	r[0].a[ta].x.j+=vec.j;
        	r[0].a[ta].x.k+=vec.k;
        	}
	return;
	}
if((xs>r[0].a[0].nalt-1)||(xd>=r[0].a[0].nalt-1)){
	mywhine("Unallocated source or destination coordinates in translate_residue_by_XYZ.");}
if(xs==-1){
	for(ta=0;ta<r[0].na;ta++){ 
        	r[0].a[ta].xa[xd].i=r[0].a[ta].x.i+vec.i;
        	r[0].a[ta].xa[xd].j=r[0].a[ta].x.j+vec.j;
        	r[0].a[ta].xa[xd].k=r[0].a[ta].x.k+vec.k;
        	}
	return;
	}
if(xd==-1){
	for(ta=0;ta<r[0].na;ta++){ 
        	r[0].a[ta].x.i=r[0].a[ta].xa[xs].i+vec.i;
        	r[0].a[ta].x.j=r[0].a[ta].xa[xs].j+vec.j;
        	r[0].a[ta].x.k=r[0].a[ta].xa[xs].k+vec.k;
        	}
	return;
	}
// still here?  
if((xs<0)||(xd<0)){mywhine("Something is very wrong in translate_residue_by_XYZ.");}
for(ta=0;ta<r[0].na;ta++){ 
       	r[0].a[ta].xa[xd].i=r[0].a[ta].xa[xs].i+vec.i;
       	r[0].a[ta].xa[xd].j=r[0].a[ta].xa[xs].j+vec.j;
       	r[0].a[ta].xa[xd].k=r[0].a[ta].xa[xs].k+vec.k;
       	}

return;
}
void translate_molecule_by_XYZ(molecule *m,int xs,int xd,coord_3D vec){
int ta=0,tr=0;
if(m[0].nr<=0){mywhine("Molecule contains no residues in translate_molecule_by_XYZ.");}
if(m[0].r[0].na<=0){mywhine("Molecule's residue zero contains no atoms in translate_molecule_by_XYZ.");}
if((xs<-1)||(xd<-1)){mywhine("Invalid source or destination coordinates to translate_molecule_by_XYZ.");}
if((xs==-1)&&(xd==-1)){
	for(tr=0;tr<m[0].nr;tr++){ 
	for(ta=0;ta<m[0].r[tr].na;ta++){ 
        	m[0].r[tr].a[ta].x.i+=vec.i;
        	m[0].r[tr].a[ta].x.j+=vec.j;
        	m[0].r[tr].a[ta].x.k+=vec.k;
        	}}
	return;
	}
if((xs>m[0].r[0].a[0].nalt-1)||(xd>=m[0].r[0].a[0].nalt-1)){
	mywhine("Unallocated source or destination coordinates in translate_molecule_by_XYZ.");}
if(xs==-1){
	for(tr=0;tr<m[0].nr;tr++){ 
	for(ta=0;ta<m[0].r[tr].na;ta++){ 
        	m[0].r[tr].a[ta].xa[xd].i=m[0].r[tr].a[ta].x.i+vec.i;
        	m[0].r[tr].a[ta].xa[xd].j=m[0].r[tr].a[ta].x.j+vec.j;
        	m[0].r[tr].a[ta].xa[xd].k=m[0].r[tr].a[ta].x.k+vec.k;
        	}}
	return;
	}
if(xd==-1){
	for(tr=0;tr<m[0].nr;tr++){ 
	for(ta=0;ta<m[0].r[tr].na;ta++){ 
        	m[0].r[tr].a[ta].x.i=m[0].r[tr].a[ta].xa[xs].i+vec.i;
        	m[0].r[tr].a[ta].x.j=m[0].r[tr].a[ta].xa[xs].j+vec.j;
        	m[0].r[tr].a[ta].x.k=m[0].r[tr].a[ta].xa[xs].k+vec.k;
        	}}
	return;
	}
// still here?  
if((xs<0)||(xd<0)){mywhine("Something is very wrong in translate_molecule_by_XYZ.");}
for(tr=0;tr<m[0].nr;tr++){ 
for(ta=0;ta<m[0].r[tr].na;ta++){ 
       	m[0].r[tr].a[ta].xa[xd].i=m[0].r[tr].a[ta].xa[xs].i+vec.i;
       	m[0].r[tr].a[ta].xa[xd].j=m[0].r[tr].a[ta].xa[xs].j+vec.j;
       	m[0].r[tr].a[ta].xa[xd].k=m[0].r[tr].a[ta].xa[xs].k+vec.k;
       	}}

return;
}
void translate_ensemble_by_XYZ(ensemble *e,int xs,int xd,coord_3D vec){
int ta=0,tr=0,tm=0;
if(e[0].nm<=0){mywhine("Ensemble contains no molecules in translate_ensemble_by_XYZ.");}
if(e[0].m[0].nr<=0){mywhine("Ensemble's molecule zero contains no residues in translate_ensemble_by_XYZ.");}
if(e[0].m[0].r[0].na<=0){mywhine("Ensemble's molecule/residue zero contains no atoms in translate_ensemble_by_XYZ.");}
if((xs<-1)||(xd<-1)){mywhine("Invalid source or destination coordinates to translate_ensemble_by_XYZ.");}
if((xs==-1)&&(xd==-1)){
	for(tm=0;tm<e[0].nm;tm++){ 
	for(tr=0;tr<e[0].m[tm].nr;tr++){ 
	for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){ 
        	e[0].m[tm].r[tr].a[ta].x.i+=vec.i;
        	e[0].m[tm].r[tr].a[ta].x.j+=vec.j;
        	e[0].m[tm].r[tr].a[ta].x.k+=vec.k;
        	}}}
	return;
	}
if((xs>e[0].m[0].r[0].a[0].nalt-1)||(xd>=e[0].m[0].r[0].a[0].nalt-1)){
	mywhine("Unallocated source or destination coordinates in translate_ensemble_by_XYZ.");}
if(xs==-1){
	for(tm=0;tm<e[0].nm;tm++){ 
	for(tr=0;tr<e[0].m[tm].nr;tr++){ 
	for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){ 
        	e[0].m[tm].r[tr].a[ta].xa[xd].i=e[0].m[tm].r[tr].a[ta].x.i+vec.i;
        	e[0].m[tm].r[tr].a[ta].xa[xd].j=e[0].m[tm].r[tr].a[ta].x.j+vec.j;
        	e[0].m[tm].r[tr].a[ta].xa[xd].k=e[0].m[tm].r[tr].a[ta].x.k+vec.k;
        	}}}
	return;
	}
if(xd==-1){
	for(tm=0;tm<e[0].nm;tm++){ 
	for(tr=0;tr<e[0].m[tm].nr;tr++){ 
	for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){ 
        	e[0].m[tm].r[tr].a[ta].x.i=e[0].m[tm].r[tr].a[ta].xa[xs].i+vec.i;
        	e[0].m[tm].r[tr].a[ta].x.j=e[0].m[tm].r[tr].a[ta].xa[xs].j+vec.j;
        	e[0].m[tm].r[tr].a[ta].x.k=e[0].m[tm].r[tr].a[ta].xa[xs].k+vec.k;
        	}}}
	return;
	}
// still here?  
if((xs<0)||(xd<0)){mywhine("Something is very wrong in translate_ensemble_by_XYZ.");}
for(tm=0;tm<e[0].nm;tm++){ 
for(tr=0;tr<e[0].m[tm].nr;tr++){ 
for(ta=0;ta<e[0].m[tm].r[tr].na;ta++){ 
       	e[0].m[tm].r[tr].a[ta].xa[xd].i=e[0].m[tm].r[tr].a[ta].xa[xs].i+vec.i;
       	e[0].m[tm].r[tr].a[ta].xa[xd].j=e[0].m[tm].r[tr].a[ta].xa[xs].j+vec.j;
       	e[0].m[tm].r[tr].a[ta].xa[xd].k=e[0].m[tm].r[tr].a[ta].xa[xs].k+vec.k;
       	}}}

return;
}
