#include "../inc/mylib.h"
#include "../inc/molecules.h"

/* 20100223 BLFoley:  When I got here, I found that these functions all 
seemed to use a variant of:

	if(foo != NULL || foo != 0x0){free(foo);}

...but I think that should be && not ||.  So, I'm gonna change them all to
be that way and if everything seems to break, I'll undo it.


sdjfdskjfh


I say this because if a system has some weird notion of what "NULL" means,
possibly because of something a programmer did, then a pointer could actually
be null, but the compiler would still try to free it, because it would never
check the second case.  And if 0x0 is always the same as NULL, then why check
both?  I'm also about to expand these considerably to reflect the new
structure design.

*/

/********** structures from parameter_sets.h ****************/

void deallocateChiralityDescription(chirality_description *cd){
 int i;
	if(cd[0].ELGEOM != NULL && cd[0].ELGEOM != 0x0){free(cd[0].ELGEOM);} 
	if(cd[0].SPGEOM != NULL && cd[0].SPGEOM != 0x0){free(cd[0].SPGEOM);} 
	for(i=0;i<cd[0].niso;i++){
		if(cd[0].iso[i] != NULL && cd[0].iso[i] != 0x0){free(cd[0].iso[i]);} 
		}
 return ;
}

void deallocateBondType(bond_type *btp){
 //Free up the description
 if(btp[0].desc != NULL && btp[0].desc != 0x0){free(btp[0].desc);}
 //Free up the first and second atom types in bond
 if(btp[0].NT[0] != NULL && btp[0].NT[0] != 0x0){free(btp[0].NT[0]);}
 if(btp[0].NT[1] != NULL && btp[0].NT[1] != 0x0){free(btp[0].NT[1]);}
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
// free up the convenience indices & pointers
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
 //Free n terms
 for(i=0;i<atp[0].nR;i++){ if(atp[0].RD[i] != NULL && atp[0].RD[i] != 0x0){free(atp[0].RD[i]);} } 
 //Print a warning if the void pointer doesn't seem to be empty
 if(atp[0].VP!= NULL && atp[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating atom_type structure with non-NULL void pointer!\n");}
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
	for(i=0;i<tps[0].na;i++) deallocateAtype(&tps[0].a[i]);
	for(i=0;i<tps[0].nr;i++) deallocateRtype(&tps[0].r[i]);
	for(i=0;i<tps[0].nm;i++) deallocateMtype(&tps[0].m[i]);
 //Print a warning if the void pointer doesn't seem to be emp if(tps[0].VP!= NULL && tps[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating types structure with non-NULL void pointer!\n");}
 return ;
}

void deallocateParameterSet(parameter_set *ps){
 int i;
 //deallocate substructures
	for(i=0;i<ps[0].nAT;i++) deallocateAtype(&ps[0].AT[i]);
	for(i=0;i<ps[0].nRT;i++) deallocateRtype(&ps[0].RT[i]);
	for(i=0;i<ps[0].nMT;i++) deallocateMtype(&ps[0].MT[i]); 
	for(i=0;i<ps[0].nBT;i++) deallocateBondType(&ps[0].BT[i]);
	for(i=0;i<ps[0].nHBT;i++) deallocateBondType(&ps[0].HBT[i]);
	for(i=0;i<ps[0].nNBT;i++) deallocateBondType(&ps[0].NBT[i]);
	for(i=0;i<ps[0].nANT;i++) deallocateAngleType(&ps[0].ANT[i]);
	for(i=0;i<ps[0].nTRT;i++) deallocateTorsionType(&ps[0].TRT[i]);
 return ;
}

/********** structures from molecules.h ****************/

void deallocateRingEnsindex(ring_ensindex *re){
 if(re[0].P != NULL && re[0].P != 0x0){free(re[0].P);} 
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
 for(i=0;i<ms[0].nmN;i++){free ms[0].mN[i]}
 for(i=0;i<ms[0].nrN;i++){free ms[0].rN[i]}
 for(i=0;i<ms[0].naN;i++){free ms[0].aN[i]}
}

void deallocateBond(bond *bnd){
 if(bnd[0].D != NULL && bnd[0].D != 0x0){free(bnd[0].D);}
 if(bnd[0].typ != NULL && bnd[0].typ != 0x0){deallocateBondType(bnd[0].typ);} 
 return ;
} 
void deallocateBondset(bondset *bst){
 int i; 
 if(bst[0].b != NULL && bst[0].b != 0x0){
  for( i = 0; i < bst[0].n; i++) deallocateBond(&bst[0].b[i]);
  free(bst[0].b);
 } 
 return ; 
}

void deallocatefoo(foo *bar){
 int i;
	if( != NULL &&  != 0x0){free();} 
	for(i=0;i< ;i++){
		if( != NULL &&  != 0x0){free();} 
		}
 //Print a warning if the void pointer doesn't seem to be empty
 if(foo[0].VP!= NULL && foo[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating foo structure with non-NULL void pointer!\n");}
 return ;
}
typedef struct {
	char isorigin; ///< Y/y/N/n -- is this the first atom in the tree?
	int ni; ///< number incoming bonds
	ensindex *i; ///< atoms making incoming bonds
	int ii; ///< index of reference incoming (outgoing if origin) bond (0 to ni-1, default=0)
	int no; ///< number outgoing bonds
	int ic; ///< index of outgoing atom used to use to set chirality
	ensindex *o; ///< atoms (no of them) making outgoing bonds
	bonded_position_set *op; ///< position information for the no outgoing bonds
	ensindex tref; ///< index of grandparent incoming atom to use to set torsion
	double tors; ///< torsion to set with respect to indexed grandparent (and ii and ic)
	int nOV; ///< number of open valences
	bonded_position_set *OV; ///< coordinate info for open valences
	int nEC; ///< number of extra coordinates (e.g., lone pairs) defined
	bonded_position_set *EC; ///< coordinate info for the extra locations
} connection_tree; ///< for non-redundant bonding descriptions. See documentation for information.

typedef struct {
	molindex s; ///< source atom in the bond
	molindex t; ///< target atom in the bond
	double o; ///< order 
	bond_type *typ; ///< pointer to type of bond
	int i; ///< index -- for example as alternative to *typ
	char *D; ///< free-form description
} molbond; ///< a bond within a molecule
typedef struct {
	int n;  // number of these
	molbond *b; // n of them
	char *D; // free-form description
} molbondset; // set of consecutive molbonds

typedef struct {
	molindex a,b,c; // three atoms in the angle
	double ang; // the angle
	angle_type *typ; // pointer to type of angle
	char *D;// free-form descriptor
	int i; // index -- for example as alternative to *typ
} angle_index; // index to atoms in an angle
typedef struct {
	molindex a,b,c,d; // four atoms defining the torsion
	double ang; // the angle
	torsion_type *typ; // pointer to type of torsion
	char *D;// free-form descriptor
	int i; // index -- for example as alternative to *typ
} torsion_index; // index to atoms in an angle

/********** structure atom *************/
typedef struct {
	int n; ///< atom number or other identifying index
	char *N; ///< atom name 
	char *T; ///< atom type name
	char *D; ///< free-form description
	double m; ///< atom mass -- units defined per program needs
	int t; ///< type number -- must correspond to assignments of "atype" (see)
	atype *typ; ///< pointer to atype structure
	molindex moli; ///< the molecule index for this atom (address)
	int nb; ///< number of bonds within this residue (deprecated, soon to go away)
	bond *b; ///< bond structures (nb of these)
	int nmb; ///< number of bonds to other residues or molecules
	molbond *mb; ///< nmb of these (pointers to the structures)
	coord_3D x; ///< atom's main coordinates 
	coord_3D xv; ///< atom's main velocity
	int nalt; ///< number of alternate coordinate sets
	coord_3D *xa; ///< nalt of alternate coords
	int nxva; ///< number of alternate velocities defined
	coord_3D *xva; ///< nxva alternate velocities
	int nvec; ///< number of vector sets
	vectormag_3D *v; ///< vector sets
	int nch; ///< number of charges
	double *ch; ///< nch charges for this atom
	int ni; ///< number of other indices
	int *i; ///< other indices, as needed (ni of these)
	int nd; ///< number of double-precision parameters
	double *d; ///< other parameters, as needed (nd of these) 
	int nensi; ///< number of ensemble indices
	ensindex *ensi; ///< list of ensemble indices
	char *sres; ///< name of some other/original residue to which this atom belonged
	int nOD; ///< number of other descriptors
	char **OD; ///< the nOD descriptors
	int nVP; ///< number of void structures
	char cID; ///< chain identifier
	void *VP; ///< nVP of these, whatever they may be
} atom; ///< an actual atom in a residue/molecule

/********** structure residue *************/
typedef struct {
	int n; ///< residue number given in input file
	char *N; ///< residue name 
	char *D; ///< free-form description for residue
	int t; ///< type number
	rtype *typ; ///< pointer to rtype structure
	molindex moli; ///< the molecule index for this residue (address -- set .a to -1 or 0) 
	int na; ///< number of atoms in residue
	double m; ///< molecular weight
	coord_3D COM; ///< center of mass for residue
	atom *a; ///< atom structures (na of these)
	connection_tree *T; ///< na of these
	int nbs; ///< number of bond sets 
	bondset *bs; ///< (consecutive bonds, use these for plotting, etc.)
	int nr; ///< number of simple rings (no cage structures, etc.)
	int nrbs; ///< number of ring bondsets defined
	bondset *rbs; ///< bondsets for rings
	int nrc; // number of ring/reference coordinate sets defined
	coord_3D *rc; ///< coordinates for ring/reference centers
	int nrp; ///< number of ring planes defined
	plane *rp; ///< equations for average/approximate/exact/etc. ring planes (where useful)
	// need something for more complex structures eventually
	int ni; ///< number of other indices
	int *i; ///< other indices, as needed (ni of these)
	int nd; ///< number of double-precision parameters
	double *d; ///< other parameters, as needed (nd of these)
	int nensi; ///< number of ensemble indices
	ensindex *ensi; ///< list of ensemble indices
	int nOD; ///< number of other descriptors
	char **OD; ///< the nOD descriptors
	int nVP; ///< number of void structures
	void *VP; ///< nVP of these, whatever they may be
} residue; ///< an actual residue in a molecule

/********** structure molecule *************/
typedef struct {
	int n; ///< molecule number given in an input file
	int i; ///< index
	char *N; ///< free-form name for molecule
	char *D; ///< free-form description for residue
	double m; ///< molecular weight
	int mi; ///< number corresponding to the index in the parent ensemble
	int Ei; ///< number corresponding to parent ensemble
	int t; ///< type number
	mtype *typ; ///< pointer to atype structure
	coord_3D COM; ///< center of mass for molecule
	int na; ///< total number of atoms in molecule
	atom **a; ///< na of these, but should point into residues
	connection_tree *aT; ///< nat of these
	int nr; ///< number of residues
        residue *r; ///< pointers to residues
	connection_tree *rT; ///< nr of these residue-level connection trees
	int nrb; ///< number of bonds between residues
	molbond *rb; ///< nrb of these descriptions of bonds 
	int nrbs; ///< number of sets of bonds between residues (for example, linear chains)
	molbondset *rbs; ///< nrbs of these sets
	int nrc; ///< number of additional reference (ring centers, for example)
	coord_3D *rc; ///< nrc of these
	int nBOX; ///< Number of box_info structures defined
	boxinfo *BOX; ///< The nBOX structures
	//coord_3D boxl,boxh; // box dimensions
	int noi; ///< number of other indices
	int *oi; ///< other indices, as needed (ni of these)
	int nd; ///< number of double-precision parameters
	double *d; ///< other parameters, as needed (nd of these)
	int nensi; ///< number of ensemble indices
	ensindex *ensi; ///< list of ensemble indices
	int nOD; ///< number of other descriptors
	char **OD; ///< the nOD descriptors
	int nVP; ///< number of void structures
	void *VP; ///< nVP of these, whatever they may be
} molecule;

/********** structure dockinfo *************/
typedef struct{
	int i; ///< index
	int n; ///< number of docked structures represented
	coord_3D RC, *TR; ///< reference coordinate ; translation for docked structure
        vectormag_3D *Q, *QN, *QN0; ///< Quaternion x,y,z,w ; Quaternion nx,ny,nz,angle
	///< NOTE!! do -not- use typical vectormag_3D functions on these structures
	int nTors; ///< Number of Torsions 
	double *Tors; ///< nTors torsions
	double *eFEB; ///< Estimated Free Energy of Binding, kcal/mol [=(1)+(3)]
	double *eKi, *Tmp; ///< Estimated Inhibition Constant, Ki & temperature, Kelvin
	double *fDE; ///< Final Docked Energy, kcal/mol [=(1)+(2)] or [=(1)+(2)+(3)-(4)]
	double *fIE; ///< (1) Final Intermolecular Energy   
	double *fIEL; ///< (2) Final Internal Energy of Ligand 
	double *TFE; ///< (3) Torsional Free Energy          
	double *USE; ///< (4) Unbound System's Energy <Autodock 4.0 only>
	int nDIH; ///< number of dihedrals
	double *DIH; ///< the dihedrals
	char *D; ///<Free-form descriptor for the set
	char *DOCK_PROGRAM; ///<The docking program used, if known and relevant
	char *VERSION; ///<The version, or similar indentifier, of the above program
	molecule M; ///< to hold molecule info about the initial/docked structures
} dockinfo;

/********** structure assembly *************/
typedef struct { 
	int i; ///< index
	char *N; ///< name
	char *D; ///< description (free-form) (was *desc)
	double mass; ///< mass of assembly
	double TIME; ///< current/initial/reference time in the simulation
	coord_3D COM; ///< center of mass 
	int nm; ///< number of molecule structures
	molecule **m; ///< nm of these
	connection_tree *mT; ///< nm of these
	//int nmb; // number of bonds/connections between molecules (H-bonds, for example)
	//molbond *mb; // nmb of these descriptions of connection 
	int nr; ///< number of residues
	residue **r; ///< nr of these -- probably best to point into molecules...
	connection_tree *rT; ///< nr of these
	int na; ///< number of atoms
	atom **a; ///< na of these -- probably best to point into mol/res combos
	int naT; ///< naT connection trees defined.
	connection_tree *aT; ///< naT of these
	int nb; ///< total number of bonds
	molbond *b; ///< the bonds
	int nANG; ///< # of angles
	angle_index *ANG; ///< nANG of these
	int nTOR; ///< # of torsions
	torsion_index *TOR; ///< nTOR of these
	int nPRM; ///< number of parameter sets
	parameter_set *PRM; ///< pointer to parameter sets
	int nmbs; ///< number of sets of connections between molecules 
	molbondset *mbs; ///< nmbs of these sets
	int nBOX; ///< Number of box_info structures defined
	boxinfo *BOX; ///< The nBOX structures
	//coord_3D boxl,boxh; // simple box dimensions (limits low & high or just use one for size)
	//double boxang; // coord_3D is x, y and z -- angle (probably between XY and YZ planes in degrees, maybe).  
	//char *boxtype; // type of simple box (cubic, trapezoidal, etc.)
	int nOD; ///< number of other descriptors
	char **OD; ///< the nOD descriptors
	int nensi; ///< number of ensemble indices
	ensindex *ensi; ///< list of ensemble indices
	int nVP; ///< number of void pointers
	void *VP; ///< void pointers
} assembly;///< structure for groups of molecules within a larger structure

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


void deallocateBondType(bond_type btp){
 //Free up the description
 if(btp.desc != NULL && btp.desc != 0x0){free(btp.desc);}
 //Free up the first and second atom types in bond
 if(btp.NT[0] != NULL && btp.NT[0] != 0x0){free(btp.NT[0]);}
 if(btp.NT[1] != NULL && btp.NT[1] != 0x0){free(btp.NT[1]);}

 return ;
}

void deallocateMolbond(molbond mlb){
 //Free up the Free Form Description
 if(mlb.D != NULL && mlb.D != 0x0){free(mlb.D);}
 //Free up the Type of Bond
 if(mlb.typ != NULL && mlb.typ != 0x0){deallocateBondType(mlb.typ[0]);}

 return ;
}

void deallocateConnectionTree(connection_tree con){
 //Free up the Incoming Bonds
// if((con.i != NULL && con.i != 0x0) && con.ni > 0){free(con.i);}
 //Free up the Outgoing Bonds
// if((con.o != NULL && con.o != 0x0) && con.no > 0){free(con.o);}
 //Free up the Description of Isomerism
// if(con.iso != NULL && con.iso != 0x0){free(con.iso);}
 //Free up the Angles
// if(con.angle != NULL && con.angle != 0x0){free(con.angle);}
 //Free up the Chirality
// if(con.chirot != NULL && con.chirot != 0x0){free(con.chirot);}
 //Free up the Distance
// if(con.distance != NULL && con.distance != 0x0){free(con.distance);}
 
 return ;
}

void deallocateBond(bond bnd){
 //Free up the Free Form Description
 if(bnd.D != NULL && bnd.D != 0x0){free(bnd.D);}
 //Free up the Type of Bond
 if(bnd.typ != NULL && bnd.typ != 0x0){deallocateBondType(bnd.typ[0]);}

 return ;
}

void deallocateBondset(bondset bst){
 int i; 
 //Free up the Bonds
 if((bst.b != NULL && bst.b != 0x0) && bst.n > 0){
  for( i = 0; i < bst.n; i++)
   deallocateBond(bst.b[i]);
  free(bst.b);
 }

 return ; 
}

void deallocateAtom(atom atm){
 int i;
 
 //Free up the Name
 if(atm.N != NULL && atm.N != 0X0){free(atm.N);}
 //Free up the Atom Type Name
 if(atm.T != NULL && atm.T != 0X0){free(atm.T);}
 //Free up the Free-Form Description
 if(atm.D != NULL && atm.D != 0X0){free(atm.D);}
 //Free up the Bond Structure
/* if((atm.b != NULL && atm.b != 0X0) && atm.nb > 0)
  for(i = 0; i < atm.nb; i++)
   deallocateMolbond(atm.b[i]);*/
 //Free up the Alternate Coordinates
 if((atm.xa != NULL && atm.xa != 0X0) && atm.nalt > 0){free(atm.xa);}
 //Free up the Vector Sets
 if((atm.v != NULL && atm.v != 0X0) && atm.nvec > 0){free(atm.v);}
 //Free up the Changes
 if((atm.ch != NULL && atm.ch != 0X0) && atm.nch > 0){free(atm.ch);}
 //Free up the Other Indicies
 if((atm.i != NULL && atm.i != 0X0) && atm.ni > 0){free(atm.i);}
 //Free up the Other Parameters
 if((atm.d != NULL && atm.d != 0X0) && atm.nd > 0){free(atm.d);}
 //Free up the list of Ensemble Indicies
 if((atm.ensi != NULL && atm.ensi != 0X0) && atm.nensi > 0){free(atm.ensi);}
 //Free up the Alternate Residue
 if(atm.sres != NULL && atm.sres != 0X0){free(atm.sres);}
 //Free up the Other Descriptors
 if((atm.OD != NULL && atm.OD != 0X0) && atm.nOD > 0){
  for(i = 0; i < atm.nOD; i++)
   free(atm.OD[i]);
  free(atm.OD);
 }
 return;
}

void deallocateResidue(residue res){
 int i;
 
 //Free up the Name
 if(res.N != NULL && res.N != 0X0){free(res.N);}
 //Free up the Free Form Description
 if(res.D != NULL && res.D != 0X0){free(res.D);}
 //Free up the Atoms
 if((res.a != NULL && res.a != 0X0) && res.na > 0){
  for(i = 0; i < res.na; i++)
   deallocateAtom(res.a[i]);
  free(res.a);
 }
 //Free up the Connection Tree
/* if((res.T != NULL && res.T != 0X0) && res.na > 0){
  for(i = 0; i < res.na; i++)
   deallocateConnectionTree(res.T[i]);
  free(res.T);
 }*/
 //Free up the Bondset
/* if((res.bs != NULL && res.bs != 0X0) && res.nbs > 0){
  for(i = 0; i < res.nbs; i++)
   deallocateBondset(res.bs[i]);
  free(res.bs);
 }*/
 //Free up the Bondset for Rings
/* if((res.rbs != NULL && res.rbs != 0X0) && res.nrbs > 0){
  for(i = 0; i < res.nrbs; i++)
   deallocateBondset(res.rbs[i]);
  free(res.rbs);
 }*/
 //Free up the Rings/Reference Centers 
 if((res.a != NULL && res.a != 0X0) && res.na > 0){free(res.rc);}
 //Free up the Ring Planes  
 if((res.rp != NULL && res.rp != 0X0) && res.nrp > 0){free(res.rp);}
 //Free up the Other Indices 
 if((res.i != NULL && res.i != 0X0) && res.ni > 0){free(res.i);}
 //Free up the Other Parameters 
 if((res.d != NULL && res.d != 0X0) && res.nd > 0){free(res.d);}
 //Free up the list of Ensemble Indices
 if((res.ensi != NULL && res.ensi != 0X0) && res.nensi > 0){free(res.ensi);}
 //Free up the Other Descriptors
 if((res.OD != NULL && res.OD != 0X0) && res.nOD > 0){
  for(i = 0; i < res.nOD; i++)
   free(res.OD[i]);
  free(res.OD);
 }
 return ;
}

void deallocateMolecule(molecule mol){
 int i;

 //Free up the Name
 if(mol.N != NULL && mol.N != 0X0){free(mol.N);}
 //Free up the Free Form Descriptor
 if(mol.D != NULL && mol.D != 0X0){free(mol.D);}
 //The atoms should be removed via their respective residues
 //Free up the Residues
 if((mol.r != NULL && mol.r != 0X0) && mol.nr > 0){
  for(i = 0; i < mol.nr; i++)
   deallocateResidue(mol.r[i]);
  free(mol.r);
 }
 //Free up the Connection Tree
/* if((mol.rT != NULL && mol.rT != 0X0) && res.nr > 0){
  for(i = 0; i < mol.nr; i++)
   deallocateConnectionTree(mol.rT[i]);
  free(mol.rT);
 }*/
 //Free up the Molbond
/* if((mol.rb != NULL && mol.rb != 0X0) && mol.nrb > 0){
  for(i = 0; i < mol.nrb; i++)
   deallocateMolbond(mol.rb[i]);
  free(mol.rb);
 }*/ 
 //Free up the Molbondset
/* if((mol.rbs != NULL && mol.rbs != 0X0) && mol.nrbs > 0){
  for(i = 0; i < mol.nrbs; i++)
   deallocateMolbondset(mol.rbs[i]);
  free(mol.rbs);
 }*/ 
 //Free up the Additional References
 if((mol.rc != NULL && mol.rc != 0X0) && mol.nrc > 0){free(mol.rc);}
 //Free up the Box Info
/* if((mol.BOX != NULL && mol.BOX != 0X0) && mol.nBOX > 0){
  for(i = 0; i < mol.nBOX; i++)
   deallocateBoxinfo(mol.Box[i]);
  free(mol.BOX);
 }*/
 //Free up the Other Indices
 if((mol.oi != NULL && mol.oi != 0X0) && mol.noi > 0){free(mol.oi);}
 //Free up the Other Parameters
 if((mol.d != NULL && mol.d != 0X0) && mol.nd > 0){free(mol.d);}
 //Free up the Ensemble Indices
 if((mol.ensi != NULL && mol.ensi != 0X0) && mol.nensi > 0){free(mol.ensi);}
 //Free up the Other Descriptors
 if((mol.OD != NULL && mol.OD != 0X0) && mol.nOD > 0){
  for(i = 0; i < mol.nOD; i++)
   free(mol.OD[i]);
  free(mol.OD);
 }
 return ;
}
