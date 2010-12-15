/** \file bonding_utilities.c  20101214 BLFoley
 * Purpose:  Provide utilities relevant to bonding within the various
 * 	structures.  For example: build a connection tree from atom-
 * 	level bonds; build atom-level bonds from a connection tree; 
 * 	deduce bonding via distance. 
 *
 * Some useful definitions:
 *
 * 	See also the molecule, residue and atom structures.  Most notation
 * 		is very consistent and can be determined from context.
 *
 * 	connection_tree:  non-redundant description of bonding
 * 	m.aT:	atom-level connectiontree in a molecule structure
 * 	m.rT:	residue-level connectiontree in a molecule structure
 * 	r.aT:	atom-level connectiontree in a residue structure
 *
 * Some caveats/requirements:
 *
 * 	* All functions here assume that the information contained in
 * 		the data they depend on are correct.
 * 	* All functions here should check for the presence of data and
 * 		write an error to stdout if the data is not found.  This
 * 		error is (hopefully) for the programmer, not end user.
 * 	* In general, all molecule indices need to be set.
 *
 * Return values:
 *
 * 	0	All is ok
 * 	1	No data present to work with (procedure state unknown)
 * 	2	Insufficient data to work with (procedure state unknown)
 * 	3	No atoms found in molecule/residue
 * 	4	No residues found in molecule
 * 
 */

/***********   follow_molecule_atom_molbonds_for_contree() ************/
/** The purpose of this function is to follow atoms recursively through
 * atom-level bonding until all are seen. 
 */
void follow_molecule_atom_molbonds_for_contree(molecule *m, moli mi){
/* If this atom has been seen in the tree, just go back */

/* Mark this atom seen in tree */

/* Announce it to all the outgoing bonds */
/* Call this function for all atoms with outgoing bonds */

return;
}

/***********   set_connection_tree_molecule_atoms() ******************/
/** The purpose of this function is to set the connection tree within a
 * molecule that is based on atom-level bonding.  
 *
 * The molecule can be alone or part of an assembly/ensemble.
 */

int set_connection_tree_molecule_atoms(molecule *m){
int 	r=0, /* residue index */
	a=0, /* atom index */
	na=0, /* number of atoms in molecule */
	*new; /* na of these; atom has NOT been seen=0, has been seen=1;

/* Do a little sanity check */
if(m[0].nr==0){
	printf("No residues found in molecule in set_connection_tree_molecule_atoms\n");
	return 4;
	}
if(m[0].r[0].na==0){
	printf("No atoms found in first residue in set_connection_tree_molecule_atoms\n");
	return 2;
	}
if(m[0].r[0].a[0].nmb==0){
	printf("No bonding info found for first atom in set_connection_tree_molecule_atoms\n");
	return 1;
	}

/* Count the number of atoms */
for(r=0;r<m[0].nr;r++){ na+=m[0].r[r].na; }
/* allocate spaces */
new=(int*)calloc(na,sizeof(int));
m[0].aT=(connection_tree*)calloc(na,sizeof(connection_tree));
/* initialize the tree */
for(a=0;a<na;a++){m[0].aT[a].ni=-1; /* indicate not seen */}
/* Start at the first atom */
m[0].aT[0].ii=-1; /* Set this as the origin atom */
m[0].aT[0].ensi.E=-1;
m[0].aT[0].ensi.A=-1;
m[0].aT[0].ensi.m=m[0].r[0].a[0].moli.m;
m[0].aT[0].ensi.r=m[0].r[0].a[0].moli.r;
m[0].aT[0].ensi.a=m[0].r[0].a[0].moli.a;

/* Call the following function starting with this one. */
follow_molecule_atom_molbonds_for_contree(m,m[0].r[0].a[0].moli);
/* That oughtta do it */

return 0;
}


/***********   set_connection_tree_molecule_residues() ******************/
/** The purpose of this function is to set the connection tree within a
 * molecule that is based on residue-level bonding.  
 */

/***************   set_connection_tree_residues_atoms() ********************/
/** The purpose of this function is to set the connection tree within a
 * residue that is based on atom-level bonding.  
 */
