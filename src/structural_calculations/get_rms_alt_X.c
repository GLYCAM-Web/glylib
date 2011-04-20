// Function written by B. Lachele Foley, starting 20071017
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* These functions determine the RMS between atoms in the 
 * relevant structures (assembly, molecule, residue, etc.).
 * This is a work in progress -- all functions may not be complete.
 *

 Functional form:

	double get_alt_rms_X(X-struct *X, int xs, int xt);

	X-struct is molecule, assembly, etc.
	X is an abbreviation (mol, res, etc.)
	xs is the location of the source coords (-1 for main, 0+ for alt)
	xt is the location of the target coords (must be an alternate coord) 
*/
/******************* get_alt_rms_res(residue *r, int xs, int xt) *****************/
double get_alt_rms_res(residue *r, int xs, int xt){
int ai=0;
double RMS=0;
coord_3D c;

if(xt<0){mywhine("get_alt_rms_res: target coord cannot be main set");}
if(xs<0){mywhine("get_alt_rms_res: main set source coords not written yet");}
for(ai=0;ai<r[0].na;ai++){
	c=subtract_coord(r[0].a[ai].xa[xs],r[0].a[ai].xa[xt]);
	RMS+=c.i*c.i+c.j*c.j+c.k*c.k;
	}
RMS/=r[0].na;
RMS=sqrt(RMS);
return RMS;
}

/******************* get_alt_rms_mol(molecule *m, int xs, int xt) *****************/
double get_alt_rms_mol(molecule *m, int xs, int xt){
int ai=0,ri=0,nat=0;
double RMS=0;
coord_3D c;

if(xt<0){mywhine("get_alt_rms_mol: target coord cannot be main set");}
if(xs<0){mywhine("get_alt_rms_mol: main set source coords not written yet");}

for(ri=0;ri<m[0].nr;ri++){
	nat+=m[0].r[ri].na;
	for(ai=0;ai<m[0].r[ri].na;ai++){
		c=subtract_coord(m[0].r[ri].a[ai].xa[xs],m[0].r[ri].a[ai].xa[xt]);
		RMS+=c.i*c.i+c.j*c.j+c.k*c.k;
	}}
RMS/=nat;
RMS=sqrt(RMS);
return RMS;
}

