/* File vector_movements.c begun on 20071231 by BLFoley.
 * Purpose: Provide functions that will move the atoms in a
 * structure according to vectors contained in the structure.
 */
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/****************** shift_molecule_atoms_by_vector() *********************/
/* Moves each atom in a molecule acording to a vector saved within the
 * molecule.  Moves the atoms according to the requested scale times the
 * vector.  For coordinates, -1 means use main coordinate set. Usage:
 * shift_molecule_atoms_by_vector(molecule *m,int xs,int xt,int vr, double scale)
 * 	molecule *m 	= the molecule containing the data
 *	int xs		= the source coordinates (the coords to be shifted)
 *	int xt		= the target coordinated (where to put the shifted coords)
 *	int vr		= the vector by which to shift the atoms
 *	double scale	= vector scale factor to use for the shifting
 */

void shift_molecule_atoms_by_vector_scale(molecule *m,int xs,int xt,int vr, double scale){
int ta=0,tb,tuse=0,tsz=0,tbins=0;
coord_3D c;

// move atoms by scale * vector
if((xs==-1)&&(xt==-1)){
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			c=scalarmult_coord(vec_to_coord(m[0].r[ta].a[tb].v[vr]),scale);
			m[0].r[ta].a[tb].x=add_coord(m[0].r[ta].a[tb].x,c);
			}
		}
	}
if((xs!=-1)&&(xt==-1)){
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			tuse=malloc_usable_size(m[0].r[ta].a[tb].xa);
			tsz=sizeof(coord_3D);
			tbins=tuse/tsz; 
			if(tbins<(xs+1)){mywhine("memory issue: source coords (1) not allocated in shift_molecule_atoms_by_vector_scale");}
			c=scalarmult_coord(vec_to_coord(m[0].r[ta].a[tb].v[vr]),scale);
			m[0].r[ta].a[tb].x=add_coord(m[0].r[ta].a[tb].xa[xs],c);
			}
		}
	}
if((xs==-1)&&(xt!=-1)){
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			tuse=malloc_usable_size(m[0].r[ta].a[tb].xa);
			tsz=sizeof(coord_3D);
			tbins=tuse/tsz; 
			if(tbins<(xt+1)){mywhine("memory issue: translated coords (1) not allocated in shift_molecule_atoms_by_vector_scale");}
			c=scalarmult_coord(vec_to_coord(m[0].r[ta].a[tb].v[vr]),scale);
			m[0].r[ta].a[tb].xa[xt]=add_coord(m[0].r[ta].a[tb].x,c);
			}
		}
	}
if((xs!=-1)&&(xt!=-1)){
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			tuse=malloc_usable_size(m[0].r[ta].a[tb].xa);
			tsz=sizeof(coord_3D); 
			tbins=tuse/tsz; 
			if(tbins<(xs+1)){mywhine("memory issue: source coords (2) not allocated in shift_molecule_atoms_by_vector_scale");}
			if(tbins<(xt+1)){mywhine("memory issue: translated coords (2) not allocated in shift_molecule_atoms_by_vector_scale");}
			c=scalarmult_coord(vec_to_coord(m[0].r[ta].a[tb].v[vr]),scale);
			m[0].r[ta].a[tb].xa[xt]=add_coord(m[0].r[ta].a[tb].xa[xs],c);
			}
		}
	}
return;
}
