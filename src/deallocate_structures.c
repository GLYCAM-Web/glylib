#include "../inc/mylib.h"
#include "../inc/molecules.h"

/* 20100223 BLFoley:  When I got here, I found that these functions all 
seemed to use a variant of:

	if(foo != NULL || foo != 0x0){free(foo);}

...but I think that should be && not ||.  So, I'm gonna change them all to
be that way and if everything seems to break, I'll undo it.

I say this because if a system has some weird notion of what "NULL"
means, possibly because of something a programmer did, then a pointer 
could actually be null, but the compiler would still try to free it,
because it would never check the second case.  And if 0x0 is always
the same as NULL, then why check both?

I'm also about to expand these considerably to reflect the new structure design.

Notes regarding these functions:

* If a structure does not contain any pointers, it can just be freed.
* Structures that contain pointers must have a deallocation function.
* Single-pointers to arrays of simple types (e.g. int) may be freed.
* Single-pointers to arrays of structures might need individual deallocation.
* Double-pointers
	* Be very careful before freeing double-pointed structures
	* 	-- they might point somewhere you don't want freed
	* Freeing the top-level pointer should not interfere with data below 
*/

/********** structures from parameter_sets.h ****************/

void deallocateChiralityDescription(chirality_description *cd){
 int i;
	if(cd[0].ELGEOM != NULL && cd[0].ELGEOM != 0x0){free(cd[0].ELGEOM);} 
	if(cd[0].SPGEOM != NULL && cd[0].SPGEOM != 0x0){free(cd[0].SPGEOM);} 
	if(cd[0].iso != NULL && cd[0].iso != 0x0){
		for(i=0;i<cd[0].niso;i++){if(cd[0].iso[i] != NULL && cd[0].iso[i] != 0x0){free(cd[0].iso[i]);}}
		free(cd[0].iso);}
 return ;
} 
void deallocateBondType(bond_type *btp){
 //Free up the description
 if(btp[0].desc != NULL && btp[0].desc != 0x0){free(btp[0].desc);}
 //Free up the first and second atom types in bond
 if(btp[0].NT[0] != NULL && btp[0].NT[0] != 0x0){free(btp[0].NT[0]);}
 if(btp[0].NT[1] != NULL && btp[0].NT[1] != 0x0){free(btp[0].NT[1]);}
 if(btp[0].NT != NULL && btp[0].NT != 0x0){free(btp[0].NT);}
 //Print a warning if the void pointer doesn't seem to be empty
 if(btp[0].VP!= NULL && btp[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating bond_type structure with non-NULL void pointer!\n");}
 return ;
}
void deallocateAngleType(angle_type *atp){
 //Free up the description
 if(atp[0].desc != NULL && atp[0].desc != 0x0){free(atp[0].desc);}
 //Free up the first and second atom types in bond
 if(atp[0].NT[0] != NULL && atp[0].NT[0] != 0x0){free(atp[0].NT[0]);}
 if(atp[0].NT[1] != NULL && atp[0].NT[1] != 0x0){free(atp[0].NT[1]);}
 if(atp[0].NT[2] != NULL && atp[0].NT[2] != 0x0){free(atp[0].NT[2]);}
 if(atp[0].NT != NULL && atp[0].NT != 0x0){free(atp[0].NT);}
 //Print a warning if the void pointer doesn't seem to be empty
 if(atp[0].VP!= NULL && atp[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating angle_type structure with non-NULL void pointer!\n");}
 return ;
}
void deallocateTorsionType(torsion_type *ttp){
 //Free up the description
 if(ttp[0].desc != NULL && ttp[0].desc != 0x0){free(ttp[0].desc);}
 //Free up the first and second atom types in bond
 if(ttp[0].NT[0] != NULL && ttp[0].NT[0] != 0x0){free(ttp[0].NT[0]);}
 if(ttp[0].NT[1] != NULL && ttp[0].NT[1] != 0x0){free(ttp[0].NT[1]);}
 if(ttp[0].NT[2] != NULL && ttp[0].NT[2] != 0x0){free(ttp[0].NT[2]);}
 if(ttp[0].NT[3] != NULL && ttp[0].NT[3] != 0x0){free(ttp[0].NT[3]);}
 if(ttp[0].NT != NULL && ttp[0].NT != 0x0){free(ttp[0].NT);}
 //Free the n terms
 if(ttp[0].k != NULL && ttp[0].k != 0x0){free(ttp[0].k);} 
 if(ttp[0].N != NULL && ttp[0].N != 0x0){free(ttp[0].N);} 
 if(ttp[0].P != NULL && ttp[0].P != 0x0){free(ttp[0].P);} 
 //Print a warning if the void pointer doesn't seem to be empty
 if(ttp[0].VP!= NULL && ttp[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating torsion_type structure with non-NULL void pointer!\n");}
 return ;
}

void deallocateAtype(atype *atp){
 int i;
 if(atp[0].N != NULL && atp[0].N != 0x0){free(atp[0].N);}
 if(atp[0].NT != NULL && atp[0].NT != 0x0){free(atp[0].NT);}
 if(atp[0].desc != NULL && atp[0].desc != 0x0){free(atp[0].desc);}
 if(atp[0].bo != NULL && atp[0].bo != 0x0){free(atp[0].bo);}
 if(atp[0].ch != NULL && atp[0].ch != 0x0){free(atp[0].ch);}
 if(atp[0].R != NULL && atp[0].R != 0x0){free(atp[0].R);}
 if(atp[0].SC != NULL && atp[0].SC != 0x0){free(atp[0].SC);}
//Free n terms
 if(atp[0].RD != NULL && atp[0].RD != 0x0){
 	for(i=0;i<atp[0].nR;i++){if(atp[0].RD[i] != NULL && atp[0].RD[i] != 0x0){free(atp[0].RD[i]);}}
	free(atp[0].RD);}
// free up the convenience indices & pointers
//	only need free top level here.  None of these should hold
//	data of their own.  They should only point into other arrays.
 if(atp[0].iBT != NULL && atp[0].iBT != 0x0){free(atp[0].iBT);}
 if(atp[0].BT != NULL && atp[0].BT != 0x0){free(atp[0].BT);}
 if(atp[0].iHBT != NULL && atp[0].iHBT != 0x0){free(atp[0].iHBT);}
 if(atp[0].HBT != NULL && atp[0].HBT != 0x0){free(atp[0].HBT);}
 if(atp[0].iNBT != NULL && atp[0].iNBT != 0x0){free(atp[0].iNBT);}
 if(atp[0].NBT != NULL && atp[0].NBT != 0x0){free(atp[0].NBT);}
 if(atp[0].iANT != NULL && atp[0].iANT != 0x0){free(atp[0].iANT);}
 if(atp[0].ANT != NULL && atp[0].ANT != 0x0){free(atp[0].ANT);}
 if(atp[0].iHANT != NULL && atp[0].iHANT != 0x0){free(atp[0].iHANT);}
 if(atp[0].HANT != NULL && atp[0].HANT != 0x0){free(atp[0].HANT);}
 if(atp[0].iNANT != NULL && atp[0].iNANT != 0x0){free(atp[0].iNANT);}
 if(atp[0].NANT != NULL && atp[0].NANT != 0x0){free(atp[0].NANT);}
 if(atp[0].iTRT != NULL && atp[0].iTRT != 0x0){free(atp[0].iTRT);}
 if(atp[0].TRT != NULL && atp[0].TRT != 0x0){free(atp[0].TRT);}
 if(atp[0].iHTRT != NULL && atp[0].iHTRT != 0x0){free(atp[0].iHTRT);}
 if(atp[0].HTRT != NULL && atp[0].HTRT != 0x0){free(atp[0].HTRT);}
 if(atp[0].iNTRT != NULL && atp[0].iNTRT != 0x0){free(atp[0].iNTRT);}
 if(atp[0].NTRT != NULL && atp[0].NTRT != 0x0){free(atp[0].NTRT);}
 //Print a warning if the void pointer doesn't seem to be empty
 if(atp[0].VP!= NULL && atp[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating atype structure with non-NULL void pointer!\n");}
 return ;
} 
void deallocateRtype(rtype *rtp){
 if(rtp[0].ac != NULL && rtp[0].ac != 0x0){free(rtp[0].ac);} 
 //Print a warning if the void pointer doesn't seem to be empty
 if(rtp[0].VP!= NULL && rtp[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating rtype structure with non-NULL void pointer!\n");}
 return ;
}
void deallocateMtype(mtype *mtp){
 //Print a warning if the void pointer doesn't seem to be empty
 if(mtp[0].VP!= NULL && mtp[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating mtype structure with non-NULL void pointer!\n");}
 return ;
} 
void deallocateTypes(types *tps){
 int i;
 //deallocate substructures
	if(tps[0].a != NULL && tps[0].a != 0x0){
		for(i=0;i<tps[0].na;i++) if(tps[0].a[i] != NULL && tps[0].a[i] != 0x0){deallocateAtype(&tps[0].a[i]);}
		free(tps[0].a);}
	if(tps[0].r != NULL && tps[0].r != 0x0){
		for(i=0;i<tps[0].nr;i++) if(tps[0].r[i] != NULL && tps[0].r[i] != 0x0){deallocateRtype(&tps[0].r[i]);}
		free(tps[0].r);}
	if(tps[0].m != NULL && tps[0].m != 0x0){
		for(i=0;i<tps[0].nm;i++) if(tps[0].m[i] != NULL && tps[0].m[i] != 0x0){deallocateMtype(&tps[0].m[i]);}
		free(tps[0].m);}
 //Print a warning if the void pointer doesn't seem to be emp if(tps[0].VP!= NULL && tps[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating types structure with non-NULL void pointer!\n");}
if(tps[0].VP!= NULL && tps[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating types structure with non-NULL void pointer!\n");}
 return ;
}

void deallocateParameterSet(parameter_set *ps){
 int i;
 //deallocate substructures
	if(ps[0].AT != NULL && ps[0].AT != 0x0){
		for(i=0;i<ps[0].nAT;i++) if(ps[0].AT[i] != NULL && ps[0].AT[i] != 0x0){deallocateAtype(&ps[0].AT[i]);}
		free(ps[0].AT);}
	if(ps[0].RT != NULL && ps[0].RT != 0x0){
		for(i=0;i<ps[0].nRT;i++) if(ps[0].RT[i] != NULL && ps[0].RT[i] != 0x0){deallocateRtype(&ps[0].RT[i]);}
		free(ps[0].RT);}
	if(ps[0].MT != NULL && ps[0].MT != 0x0){
		for(i=0;i<ps[0].nMT;i++) if(ps[0].MT[i] != NULL && ps[0].MT[i] != 0x0){deallocateMtype(&ps[0].MT[i]); }
		free(ps[0].MT);}
	if(ps[0].BT != NULL && ps[0].BT != 0x0){
		for(i=0;i<ps[0].nBT;i++) if(ps[0].BT[i] != NULL && ps[0].BT[i] != 0x0){deallocateBondType(&ps[0].BT[i]);}
		free(ps[0].BT);}
	if(ps[0].HBT != NULL && ps[0].HBT != 0x0){
		for(i=0;i<ps[0].nHBT;i++) if(ps[0].HBT[i] != NULL && ps[0].HBT[i] != 0x0){deallocateBondType(&ps[0].HBT[i]);}
		free(ps[0].HBT);}
	if(ps[0].NBT != NULL && ps[0].NBT != 0x0){
		for(i=0;i<ps[0].nNBT;i++) if(ps[0].NBT[i] != NULL && ps[0].NBT[i] != 0x0){deallocateBondType(&ps[0].NBT[i]);}
		free(ps[0].NBT);}
	if(ps[0].ANT != NULL && ps[0].ANT != 0x0){
		for(i=0;i<ps[0].nANT;i++) if(ps[0].ANT[i] != NULL && ps[0].ANT[i] != 0x0){deallocateAngleType(&ps[0].ANT[i]);}
		free(ps[0].ANT);}
	if(ps[0].TRT != NULL && ps[0].TRT != 0x0){
		for(i=0;i<ps[0].nTRT;i++) if(ps[0].TRT[i] != NULL && ps[0].TRT[i] != 0x0){deallocateTorsionType(&ps[0].TRT[i]);}
		free(ps[0].TRT);}
 return ;
}

/********** structures from molecules.h ****************/

void deallocateRingEnsindex(ring_ensindex *re){
 if(re[0].P != NULL && re[0].P != 0x0){free(re[0].P);} //ensindex is simple -- no deallocation
 if(re[0].in != NULL && re[0].in != 0x0){free(re[0].in);} 
 if(re[0].out != NULL && re[0].out != 0x0){free(re[0].out);} 
 return ;
}
void deallocateResidueTreeIndex(residue_tree_index *rti){
 if(rti[0].ai != NULL && rti[0].ai != 0x0){free(rti[0].ai);} 
}
void deallocateMoleculeTreeIndex(molecule_tree_index *mti){
 if(mti[0].ri != NULL && mti[0].ri != 0x0){free(mti[0].ri);} 
 if(mti[0].r != NULL && mti[0].r != 0x0){deallocateResidueTreeIndex(&mti[0].r);} 
}
void deallocateEnsembleTreeIndex(ensemble_tree_index *eti){
 if(eti[0].mi != NULL && eti[0].mi != 0x0){free(eti[0].mi);} 
 if(eti[0].m != NULL && eti[0].m != 0x0){deallocateResidueTreeIndex(&eti[0].m);} 
}

void deallocateMoietySelection(moiety_selection *ms){
 int i;
 if(ms[0].mn != NULL && ms[0].mn != 0x0){free(ms[0].mn);} 
 if(ms[0].mi != NULL && ms[0].mi != 0x0){free(ms[0].mi);} 
 if(ms[0].rn != NULL && ms[0].rn != 0x0){free(ms[0].rn);} 
 if(ms[0].ri != NULL && ms[0].ri != 0x0){free(ms[0].ri);} 
 if(ms[0].an != NULL && ms[0].an != 0x0){free(ms[0].an);} 
 if(ms[0].ai != NULL && ms[0].ai != 0x0){free(ms[0].ai);} 
 if(ms[0].mN != NULL && ms[0].mN != 0x0){
 	for(i=0;i<ms[0].nmN;i++){if(ms[0].mN[i] != NULL && ms[0].mN[i] != 0x0){free ms[0].mN[i]}}
	free(ms[0].mN);}
 if(ms[0].rN != NULL && ms[0].rN != 0x0){
 	for(i=0;i<ms[0].nrN;i++){if(ms[0].rN[i] != NULL && ms[0].rN[i] != 0x0){free ms[0].rN[i]}}
	free(ms[0].rN);}
 if(ms[0].aN != NULL && ms[0].aN != 0x0){
 	for(i=0;i<ms[0].naN;i++){if(ms[0].aN[i] != NULL && ms[0].aN[i] != 0x0){free ms[0].aN[i]}}
	free(ms[0].aN);}
}

void deallocateBond(bond *bnd){
 if(bnd[0].D != NULL && bnd[0].D != 0x0){free(bnd[0].D);}
 // this should just be a pointer into an array, so set to null
 bnd[0].typ = 0x0; 
 return ;
} 
void deallocateBondset(bondset *bst){
 int i; 
 if(bst[0].b != NULL && bst[0].b != 0x0){
  for(i=0;i<bst[0].n;i++){deallocateBond(&bst[0].b[i]);}
  free(bst[0].b);}
 return ; 
} 
void deallocateConnectionTree(connection_tree *ct){
// all members are simple, no deallocation
 if(ct[0].i != NULL && ct[0].i != 0x0){free(ct[0].i);}
 if(ct[0].o != NULL && ct[0].o != 0x0){free(ct[0].o);}
 if(ct[0].op != NULL && ct[0].op != 0x0){free(ct[0].op);}
 if(ct[0].OV != NULL && ct[0].OV != 0x0){free(ct[0].OV);}
 if(ct[0].EC != NULL && ct[0].EC != 0x0){free(ct[0].EC);}
 return ;
}

void deallocateMolbond(molbond *mb){
 // the bond_type pointer should point into an array, so just set null
 mb[0].typ = 0x0;
 if(mb[0].D != NULL && mb[0].D != 0x0){free(mb[0].D);}
 return ;
}
void deallocateMolbondset(molbondset *mbs){
 int i;
 if(mbs[0].D != NULL && mbs[0].D != 0x0){free(mbs[0].D);}
 if(mbs[0].b != NULL && mbs[0].b != 0x0){
 	for(i=0;i<mbs[0].n;i++){deallocateMolbond(&mbs[0].b[i]);}
 	free(mbs[0].b);}
 return ;
}
void deallocateMolbondset(molbondset *mbs){
 int i;
 if(mbs[0].b != NULL && mbs[0].b != 0x0){
 	for(i=0;i<mbs[0].n;i++){deallocateMolbond(&mbs[0].b[i]);}
 	free(mbs[0].b);}
 if(mbs[0].D != NULL && mbs[0].D != 0x0){free(mbs[0].D);}
 return ;
}
void deallocateAngleIndex(angle_index *ai){
 // the pointer should point into an array, so just set null
 ai[0].typ = 0x0;
 if(ai[0].D != NULL && ai[0].D != 0x0){free(ai[0].D);}
 return ;
}
void deallocateAngleIndex(angle_index *ai){
 // the pointer should point into an array, so just set null
 ai[0].typ = 0x0;
 if(ai[0].D != NULL && ai[0].D != 0x0){free(ai[0].D);}
 return ;
}
void deallocateTorsinIndex(torsion_index *ai){
 // the pointer should point into an array, so just set null
 ti[0].typ = 0x0;
 if(ti[0].D != NULL && ti[0].D != 0x0){free(ti[0].D);}
 return ;
}

/********** structure atom *************/
void deallocateAtom(atom *a){
 int i;
 //set null
 a[0].typ = 0x0;
 // deallocate each
 if(a[0].b != NULL && a[0].b != 0x0){
 	for(i=0;i<a[0].nb;i++){deallocateBond(&a[0].b[i]);}
 	free(a[0].b);}
 if(a[0].mb != NULL && a[0].mb != 0x0){
 	for(i=0;i<a[0].nmb;i++){deallocateMolbond(&a[0].mb[i]);}
 	free(a[0].mb);}
 // free each
 if(a[0].OD != NULL && a[0].OD != 0x0){
 	for(i=0;i<a[0].nOD;i++){if(a[0].OD[i] != NULL && a[0].OD[i] != 0x0){free a[0].OD[i]}}
	free(a[0].OD);}
 // free top
 if(a[0].N != NULL && a[0].N != 0x0){free(a[0].N);}
 if(a[0].T != NULL && a[0].T != 0x0){free(a[0].T);}
 if(a[0].D != NULL && a[0].D != 0x0){free(a[0].D);}
 if(a[0].xa != NULL && a[0].xa != 0x0){free(a[0].xa);}
 if(a[0].xva != NULL && a[0].xva != 0x0){free(a[0].xva);}
 if(a[0].v != NULL && a[0].v != 0x0){free(a[0].v);}
 if(a[0].ch != NULL && a[0].ch != 0x0){free(a[0].ch);}
 if(a[0].i != NULL && a[0].i != 0x0){free(a[0].i);}
 if(a[0].d != NULL && a[0].d != 0x0){free(a[0].d);}
 if(a[0].ensi != NULL && a[0].ensi != 0x0){free(a[0].ensi);}
 if(a[0].sres != NULL && a[0].sres != 0x0){free(a[0].sres);}
	// check and warn if non-null
 if(a[0].VP!= NULL && a[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating atom structure with non-NULL void pointer!\n");}
 return ;
}
void deallocateResidue(residue *r){
 int i;
 //set null
 r[0].typ = 0x0;
 // deallocate each
 if(r[0].a != NULL && r[0].a != 0x0){
 	for(i=0;i<r[0].na;i++){deallocateAtom(&r[0].a[i]);}
 	free(r[0].a);}
 if(r[0].a != NULL && r[0].a != 0x0){
 	for(i=0;i<r[0].na;i++){deallocateAtom(&r[0].a[i]);}
 	free(r[0].a);}
 if(r[0].aT != NULL && r[0].aT != 0x0){
 	for(i=0;i<r[0].na;i++){deallocateConnectionTree(&r[0].aT[i]);}
 	free(r[0].aT);}
 if(r[0].bs != NULL && r[0].bs != 0x0){
 	for(i=0;i<r[0].nbs;i++){deallocateBondset(&r[0].bs[i]);}
 	free(r[0].bs);}
 if(r[0].rbs != NULL && r[0].rbs != 0x0){
 	for(i=0;i<r[0].nrbs;i++){deallocateBondset(&r[0].rbs[i]);}
 	free(r[0].rbs);}
 // free each
 if(r[0].OD != NULL && r[0].OD != 0x0){
 	for(i=0;i<r[0].nOD;i++){if(r[0].OD[i] != NULL && r[0].OD[i] != 0x0){free r[0].OD[i]}}
	free(r[0].OD);}
 // free top
 if(r[0].N != NULL && r[0].N != 0x0){free(r[0].N);}
 if(r[0].T != NULL && r[0].T != 0x0){free(r[0].T);}
 if(r[0].D != NULL && r[0].D != 0x0){free(r[0].D);}
 if(r[0].rc != NULL && r[0].rc != 0x0){free(r[0].rc);}
 if(r[0].rp != NULL && r[0].rp != 0x0){free(r[0].rp);}
 if(r[0].i != NULL && r[0].i != 0x0){free(r[0].i);}
 if(r[0].d != NULL && r[0].d != 0x0){free(r[0].d);}
 if(r[0].ensi != NULL && r[0].ensi != 0x0){free(r[0].ensi);}
	// check and warn if non-null
 if(r[0].VP!= NULL && r[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating residue structure with non-NULL void pointer!\n");}
 return ;
}
void deallocateMolecule(molecule *m){
 int i;
 //set null
 m[0].typ = 0x0;
 // deallocate each
 if(m[0].r != NULL && m[0].r != 0x0){
 	for(i=0;i<m[0].nr;i++){deallocateResidue(&m[0].r[i]);}
 	free(m[0].r);}
 if(m[0].rb != NULL && m[0].rb != 0x0){
 	for(i=0;i<m[0].nrb;i++){deallocateMolbond(&m[0].rb[i]);}
 	free(m[0].rb);}
 if(m[0].rbs != NULL && m[0].rbs != 0x0){
 	for(i=0;i<m[0].nrbs;i++){deallocateMolbondset(&m[0].rbs[i]);}
 	free(m[0].rbs);}
 if(m[0].aT != NULL && m[0].aT != 0x0){
 	for(i=0;i<m[0].na;i++){deallocateConnectionTree(&m[0].aT[i]);}
 	free(m[0].aT);}
 if(m[0].rT != NULL && m[0].rT != 0x0){
 	for(i=0;i<m[0].nr;i++){deallocateConnectionTree(&m[0].rT[i]);}
 	free(m[0].rT);}
 if(m[0].BOX != NULL && m[0].BOX != 0x0){
 	for(i=0;i<m[0].nBOX;i++){deallocateBoxinfo(&m[0].BOX[i]);}
 	free(m[0].BOX);}
 // free each
 if(m[0].OD != NULL && m[0].OD != 0x0){
 	for(i=0;i<m[0].nOD;i++){if(m[0].OD[i] != NULL && m[0].OD[i] != 0x0){free m[0].OD[i]}}
	free(m[0].OD);}
 // free top
 if(m[0].N != NULL && m[0].N != 0x0){free(m[0].N);}
 if(m[0].D != NULL && m[0].D != 0x0){free(m[0].D);}
 if(m[0].T != NULL && m[0].T != 0x0){free(m[0].T);}
 if(m[0].a != NULL && m[0].a != 0x0){free(m[0].a);}
 if(m[0].rc != NULL && m[0].rc != 0x0){free(m[0].rc);}
 if(m[0].oi != NULL && m[0].oi != 0x0){free(m[0].oi);}
 if(m[0].d != NULL && m[0].d != 0x0){free(m[0].d);}
 if(m[0].ensi != NULL && m[0].ensi != 0x0){free(m[0].ensi);}
	// check and warn if non-null
 if(m[0].VP!= NULL && m[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating molecule structure with non-NULL void pointer!\n");}
 return ;
} 
void deallocateDockinfo(dockinfo *di){
 int i;
 // deallocate each
 deallocateMolecule(&di[0].M); // this can't be null...
 // free top
 if(di[0].TR != NULL && di[0].TR != 0x0){free(di[0].TR);}
 if(di[0].Q != NULL && di[0].Q != 0x0){free(di[0].Q);}
 if(di[0].QN != NULL && di[0].QN != 0x0){free(di[0].QN);}
 if(di[0].QN0 != NULL && di[0].QN0 != 0x0){free(di[0].QN0);}
 if(di[0].Tors != NULL && di[0].Tors != 0x0){free(di[0].Tors);}
 if(di[0].eFEB != NULL && di[0].eFEB != 0x0){free(di[0].eFEB);}
 if(di[0].eKi != NULL && di[0].eKi != 0x0){free(di[0].eKi);}
 if(di[0].Tmp != NULL && di[0].Tmp != 0x0){free(di[0].Tmp);}
 if(di[0].fDE != NULL && di[0].fDE != 0x0){free(di[0].fDE);}
 if(di[0].fIE != NULL && di[0].fIE != 0x0){free(di[0].fIE);}
 if(di[0].fIEL != NULL && di[0].fIEL != 0x0){free(di[0].fIEL);}
 if(di[0].TFE != NULL && di[0].TFE != 0x0){free(di[0].TFE);}
 if(di[0].USE != NULL && di[0].USE != 0x0){free(di[0].USE);}
 if(di[0].DIH != NULL && di[0].DIH != 0x0){free(di[0].DIH);}
 if(di[0].D != NULL && di[0].D != 0x0){free(di[0].D);}
 if(di[0].DOCK_PROGRAM != NULL && di[0].DOCK_PROGRAM != 0x0){free(di[0].DOCK_PROGRAM);}
 if(di[0].VERSION != NULL && di[0].VERSION != 0x0){free(di[0].VERSION);}
 return ;
}
 
void deallocateAssembly(assembly *a){
 int i;
 // deallocate each
 if(a[0].mT != NULL && a[0].mT != 0x0){
 	for(i=0;i<a[0].nm;i++){deallocateConnectionTree(&a[0].mT[i]);}
 	free(a[0].mT);}
 if(a[0].RT != NULL && a[0].RT != 0x0){
 	for(i=0;i<a[0].nr;i++){deallocateConnectionTree(&a[0].RT[i]);}
 	free(a[0].RT);}
 if(a[0].aT != NULL && a[0].aT != 0x0){
 	for(i=0;i<a[0].na;i++){deallocateConnectionTree(&a[0].aT[i]);}
 	free(a[0].aT);}
 if(a[0].b != NULL && a[0].b != 0x0){
 	for(i=0;i<a[0].nb;i++){deallocateMolbond(&a[0].b[i]);}
 	free(a[0].b);}
 if(a[0].bs != NULL && a[0].bs != 0x0){
 	for(i=0;i<a[0].nbs;i++){deallocateMolbondset(&a[0].bs[i]);}
 	free(a[0].bs);}
 if(a[0].ANG != NULL && a[0].ANG != 0x0){
 	for(i=0;i<a[0].nANG;i++){deallocateAngleIndex(&a[0].ANG[i]);}
 	free(a[0].ANG);}
 if(a[0].TOR != NULL && a[0].TOR != 0x0){
 	for(i=0;i<a[0].nTOR;i++){deallocateTorsionIndex(&a[0].TOR[i]);}
 	free(a[0].TOR);}
 if(a[0].BOX != NULL && a[0].BOX != 0x0){
 	for(i=0;i<a[0].nBOX;i++){deallocateBoxinfo(&a[0].BOX[i]);}
 	free(a[0].BOX);}
 if(a[0].PRM != NULL && a[0].PRM != 0x0){
 	for(i=0;i<a[0].nPRM;i++){deallocateParameterSet(&a[0].PRM[i]);}
 	free(a[0].PRM);}
 // free each
 if(a[0].OD != NULL && a[0].OD != 0x0){
 	for(i=0;i<a[0].nOD;i++){if(a[0].OD[i] != NULL && a[0].OD[i] != 0x0){free a[0].OD[i]}}
	free(a[0].OD);}
 // free top
 if(a[0].N != NULL && a[0].N != 0x0){free(a[0].N);}
 if(a[0].T != NULL && a[0].T != 0x0){free(a[0].T);}
 if(a[0].D != NULL && a[0].D != 0x0){free(a[0].D);}
 if(a[0].r != NULL && a[0].r != 0x0){free(a[0].r);}
 if(a[0].a != NULL && a[0].a != 0x0){free(a[0].a);}
 if(a[0].m != NULL && a[0].m != 0x0){free(a[0].m);}
 if(a[0].ensi != NULL && a[0].ensi != 0x0){free(a[0].ensi);}
	// check and warn if non-null
 if(a[0].VP!= NULL && a[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating molecule structure with non-NULL void pointer!\n");}
 return ;
}
void deallocateFullAssembly(assembly *a){
 if(a[0].m != NULL && a[0].m != 0x0){
	for(i=0;i<a[0].nm;i++){if(a[0].m[i] != NULL && a[0].m[i] != 0x0){deallocateMolecule(a[0].m[i]);}}}
 deallocateAssembly(a);
 return ;
}

void deallocateResidue(residue *r){
 int i;
 //set null
 r[0].typ = 0x0;
 // deallocate each
 if(r[0].a != NULL && r[0].a != 0x0){
 	for(i=0;i<r[0].na;i++){deallocateAtom(&r[0].a[i]);}
 	free(r[0].a);}
 if(r[0].a != NULL && r[0].a != 0x0){
 	for(i=0;i<r[0].na;i++){deallocateAtom(&r[0].a[i]);}
 	free(r[0].a);}
 if(r[0].aT != NULL && r[0].aT != 0x0){
 	for(i=0;i<r[0].na;i++){deallocateConnectionTree(&r[0].aT[i]);}
 	free(r[0].aT);}
 if(r[0].bs != NULL && r[0].bs != 0x0){
 	for(i=0;i<r[0].nbs;i++){deallocateBondset(&r[0].bs[i]);}
 	free(r[0].bs);}
 if(r[0].rbs != NULL && r[0].rbs != 0x0){
 	for(i=0;i<r[0].nrbs;i++){deallocateBondset(&r[0].rbs[i]);}
 	free(r[0].rbs);}
 // free each
 if(r[0].OD != NULL && r[0].OD != 0x0){
 	for(i=0;i<r[0].nOD;i++){if(r[0].OD[i] != NULL && r[0].OD[i] != 0x0){free r[0].OD[i]}}
	free(r[0].OD);}
 // free top
 if(r[0].N != NULL && r[0].N != 0x0){free(r[0].N);}
 if(r[0].T != NULL && r[0].T != 0x0){free(r[0].T);}
 if(r[0].D != NULL && r[0].D != 0x0){free(r[0].D);}
 if(r[0].rc != NULL && r[0].rc != 0x0){free(r[0].rc);}
 if(r[0].rp != NULL && r[0].rp != 0x0){free(r[0].rp);}
 if(r[0].i != NULL && r[0].i != 0x0){free(r[0].i);}
	// check and warn if non-null
 if(m[0].VP!= NULL && m[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating molecule structure with non-NULL void pointer!\n");}
 return ;
}
/********** structure ensemble *************/
typedef struct { 
	int i; ///< index
	char *N; ///< name
	char *D; ///< free-form descriptor
	double mass; ///< mass of ensemble
	coord_3D COM; ///< center of mass 
	int nm; ///< number of molecule structures
	molecule *m; ///< nm of these
	int nA; ///< number of assembly structures 
	assembly **A; ///< nA of these
	int nBOX; ///< Number of box_info structures defined
	boxinfo *BOX; ///< The nBOX structures
	//coord_3D boxl,boxh; // box dimensions
	int nPRM; ///< number of parameter sets
	parameter_set *PRM; ///< pointer to parameter sets
	int nOD; ///< number of other descriptors
	char **OD; ///< the nOD descriptors
	int nensi; ///< number of ensemble indices
	ensindex *ensi; ///< list of ensemble indices
	int nVP; ///< number of void pointers
	void *VP; ///< void pointers
} ensemble;///< structure for an entire system of molecules

