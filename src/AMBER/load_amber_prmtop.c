/* File load_amber_prmtop.c begun on 20080608 by BLFoley
	Purpose:  read amber-style prmtop files into an
		assembly for later use. */

#include "AMBER/amber_prmtop.h"

assembly load_amber_prmtop(fileset F){
amber_prmtop *P;
assembly A;

// initialize the amber prmtop info
P=(amber_prmtop *)calloc(1,sizeof(amber_prmtop));
amber_prmtop_init(P); 

// load info as-is into sections
read_amber_prmtop_asis(F,P);

// parse sections into a glylib assembly structure
A=parse_amber_prmtop(P);

return A;
}
