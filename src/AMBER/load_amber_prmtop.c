/* File load_amber_prmtop.c begun on 20080608 by BLFoley
	Purpose:  read amber-style prmtop files into an
		assembly for later use. */

#include "AMBER/amber_prmtop.h"

assembly load_amber_prmtop(fileset F){
amber_prmtop *P;
int m,r;
assembly A;

/* initialize the amber prmtop info */
P=(amber_prmtop *)calloc(1,sizeof(amber_prmtop));
amber_prmtop_init(P); 

/* load info as-is into sections */
read_amber_prmtop_asis(F,P);

/* parse sections into a glylib assembly structure */
A=parse_amber_prmtop(P);

/* Make sure the molecule indexes are all set */
set_assembly_molindexes(&A);
/* Set the non-redundant bonding */
for(m=0;m<A.nm;m++)
	{
	for(r=0;r<A.m[m][0].nr;r++)
		{
		set_residue_atom_nodes_from_bonds(&A.m[m][0].r[r]);
		}
	set_molecule_residue_molbonds(A.m[m]);
	set_molecule_residue_nodes_from_bonds(A.m[m]);
	set_molecule_atom_nodes_from_bonds(A.m[m]);
	}

return A;
}
