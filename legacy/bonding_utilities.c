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
void follow_molecule_atom_molbonds_for_contree(molecule *m, molindex mi){
int t=m[0].r[mi.r].a[mi.a].mTi, /* the tree index */
	ii=0, /* local molbond index of incoming atom */
	rbi=0, /* convenience holder for local rbi */
	mbi=0; /* counter for molbonds */ 
molindex fi; /* convenience molindex to follow */

/* If this atom has been seen in the tree, just go back. */
if(m[0].aT[t].rbi!=-1){return;};

/* Mark this atom seen in tree */
if(m[0].aT[t].ni>0){/* if there is an incoming bond */
	rbi=m[0].aT[t].rbi=m[0].aT[t].i[0].i; /* set the index to the first one */
	/* Identify the incoming bond */
	ii=-1;
	for(mbi=0;mbi<m[0].r[mi.r].a[mi.a].nmb;mbi+++){
		fi=m[0].r[mi.r].a[mi.a].mb[mbi].t;
		if(m[0].r[fi.r].a[fi.a].moli.m!=m[0].aT[t].ensi.m){
			mywhine("molecule indices do not match in follow_molecule_atom_molbonds_for_contree");}
		if((fi.r==m[0].aT[rbi].ensi.r)&&(fi.a==m[0].aT[rbi].ensi.a)){ii=mbi;}
		}
	if(ii==-1){mywhine("did not find bond match in follow_molecule_atom_molbonds_for_contree.");}
	} 
else{ /* if this is an origin atom */
	if(m[0].r[mi.r].a[mi.a].mb==NULL){mywhine("cannot access molbonds (mol/atom/contree)");}
	ii=0;
	fi=m[0].r[mi.r].a[mi.a].mb[0].t; /* find moli of first bond target */
	rbi=m[0].aT[t].rbi=m[0].r[fi.r].a[fi.a].mTi;
	}

/*  Determine the proper bond order to maintain chirality. 

	Get each cross product vector for that with the incoming
	bond.  Save it.  Get all other cross-products.
	Get angles.  Sort angles.  Set outgoing in order.
*/

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
int 	r=ri=0, /* residue index */
	a=ai=0, /* atom index */
	t=0, /* tree index */
	na=0, /* number of atoms in molecule */
	*new; /* na of these; atom has NOT been seen=0, has been seen=1;*/

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
t=0;
ri=ai=-1;
for(r=0;r<m[0].nr;r++){ 
	for(a=0;a<m[0].r[r].na;a++){ 
		m[0].aT[t].rbi=-1; /* indicate not seen */
		m[0],aT[t].isorigin='N'; /* indicate not origin */
		m[0].aT[t].ensi.i=m[0].r[r].a[a].mTi=t; /* know thyself */
		m[0].aT[t].ensi.E=-1;
		m[0].aT[t].ensi.A=-1;
		m[0].aT[t].ensi.m=m[0].r[r].a[a].moli.m;
		m[0].aT[t].ensi.r=m[0].r[r].a[a].moli.r;
		m[0].aT[t].ensi.a=m[0].r[r].a[a].moli.a;
		/* while we're here, find the first atom that has only one bond.
			We'll use that as the tree origin */
		if(ai==-1) {
			if(m[0].r[r].a[a].nmb==1){
				m[0],aT[t].isorigin='Y'; /* indicate is origin */
				ai=a;
				ri=r;
				}
			}
		t++;
		}
	if(t>na){mywhine("tree count got larger than na in set_connection_tree_molecule_atoms");}
	}
if(na==1){ /* just a single atom */
	m[0].aT[0].ni=m[0].aT[0].no=m[0].aT[0].rbi=0;
	return 0;
	}
if(ai==-1){ /* if all atoms have more than one bond */ 
	m[0],aT[0].isorigin='Y'; /* indicate is origin */
	ai=0;
	ri=0;
	}
/* Start at the origin atom */
follow_molecule_atom_molbonds_for_contree(m,m[0].r[ri].a[ai].moli);
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
