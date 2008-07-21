/* File get_dipole.c begun on 20071231 by BLFoley
 * Purpose: contains functions for calculating dipole moments on
 * various structures.  Dipole moments should be calculated via
 * the center of mass of the unit (only matters if the unit is 
 * charged, but might as well).  Since alternate coordinate sets are
 * allowed, the COM will be calculated here (despite anything listed 
 * within the structure).
 */
#include <mylib.h>
#include <molecules.h>
/****************** get_molecule_point_charge_dipole() *********************/
/* Determines the dipole for a molecule based on the point charges of the 
 * atoms at the coordinates specified.  The atom type array is needed to 
 * supply atom masses (to get the COM for the specified coordinates).  
 */

vectormag_3D get_molecule_point_charge_dipole(molecule *m,int xs,int chgsrc,atype *AT){
int ta=0,tb,tuse=0,tsz=0,tbins=0;
coord_3D com;
vectormag_3D dip;

dip = zero_vec();
// determine the center of mass 
com=get_molecule_COM(m,AT,xs);
// calculate the dipole relative to the center of mass
if(xs==-1){
// sum over all charge*(x-COM)
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			dip.i+=m[0].r[ta].a[tb].ch[chgsrc]*(m[0].r[ta].a[tb].x.i-com.i);
			dip.j+=m[0].r[ta].a[tb].ch[chgsrc]*(m[0].r[ta].a[tb].x.j-com.j);
			dip.k+=m[0].r[ta].a[tb].ch[chgsrc]*(m[0].r[ta].a[tb].x.k-com.k);
			}
        	}
	}
else{
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			tuse=malloc_usable_size(m[0].r[ta].a[tb].xa);
                	tsz=sizeof(coord_3D);
                	tbins=tuse/tsz; 
			if(tbins<(xs+1)){mywhine("memory issue: source coords not allocated in get_molecule_point_charge_dipole");}
			dip.i+=m[0].r[ta].a[tb].ch[chgsrc]*(m[0].r[ta].a[tb].xa[xs].i-com.i);
			dip.j+=m[0].r[ta].a[tb].ch[chgsrc]*(m[0].r[ta].a[tb].xa[xs].j-com.j);
			dip.k+=m[0].r[ta].a[tb].ch[chgsrc]*(m[0].r[ta].a[tb].xa[xs].k-com.k);
			}
        	}
	}
dip.d=get_magnitude(dip);

return dip;
}
