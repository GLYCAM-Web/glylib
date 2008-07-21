// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/****************** translate_zero_to_coord_M() *********************/
/* Translates the specifiec coords (location xs) so that the coordinate
indicated by c is at the origin.  Places the translated coords in xt. */

// START HERE -- integrate or retire this ASAP

void translate_zero_to_coord_M(molecule *m,int xs,int xt,coord_3D c){
int ta=0,tb,tuse=0,tsz=0,tbins=0;

// move center of mass to specified location
if((xs==-1)&&(xt==-1)){
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			m[0].r[ta].a[tb].x=subtract_coord(m[0].r[ta].a[tb].x,c);
			}
        	}
	}
if((xs!=-1)&&(xt==-1)){
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			tuse=malloc_usable_size(m[0].r[ta].a[tb].xa);
                	tsz=sizeof(coord_3D);
                	tbins=tuse/tsz; 
			if(tbins<(xs+1)){mywhine("memory issue: source coords (1) not allocated in translate_zero_to_coord");}
			m[0].r[ta].a[tb].x=subtract_coord(m[0].r[ta].a[tb].xa[xs],c);
			}
        	}
	}
if((xs==-1)&&(xt!=-1)){
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			tuse=malloc_usable_size(m[0].r[ta].a[tb].xa);
                	tsz=sizeof(coord_3D);
                	tbins=tuse/tsz; 
			if(tbins<(xt+1)){mywhine("memory issue: translated coords (1) not allocated in translate_zero_to_coord");}
			m[0].r[ta].a[tb].xa[xt]=subtract_coord(m[0].r[ta].a[tb].x,c);
			}
        	}
	}
if((xs!=-1)&&(xt!=-1)){
	for(ta=0;ta<m[0].nr;ta++){ 
		for(tb=0;tb<m[0].r[ta].na;tb++){ 
			tuse=malloc_usable_size(m[0].r[ta].a[tb].xa);
                	tsz=sizeof(coord_3D); 
                	tbins=tuse/tsz; 
			if(tbins<(xs+1)){mywhine("memory issue: source coords (2) not allocated in translate_zero_to_coord");}
			if(tbins<(xt+1)){mywhine("memory issue: translated coords (2) not allocated in translate_zero_to_coord");}
			m[0].r[ta].a[tb].xa[xt]=subtract_coord(m[0].r[ta].a[tb].xa[xs],c);
			}
        	}
	}
return;
}
