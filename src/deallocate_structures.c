#include "../inc/mylib.h"
#include "../inc/molecules.h"

/* 20100223 BLFoley:  When I got here, I found that these functions all 
seemed to use a variant of:

	if(foo != NULL || foo != 0x0){free(foo);}

...but I think that should be && not ||.  So, I'm gonna change them all to
be that way and if everything seems to break, I'll undo it.

I'm also about to expand these considerably to reflect the new structure design.

*/

/********** structures from parameter_sets.h ****************/

void deallocateChiralityDescription(chirality_description *cd){
 int i;
	if(cd.ELGEOM != NULL && cd.ELGEOM != 0x0){free(cd.ELGEOM);} 
	if(cd.SPGEOM != NULL && cd.SPGEOM != 0x0){free(cd.SPGEOM);} 
	for(i=0;i<cd.niso;i++){
		if(cd.iso[i] != NULL && cd.iso[i] != 0x0){free(cd.iso[i]);} 
		}
 return ;
}

void deallocateBondType(bond_type *btp){
 //Free up the description
 if(btp.desc != NULL && btp.desc != 0x0){free(btp.desc);}
 //Free up the first and second atom types in bond
 if(btp.NT[0] != NULL && btp.NT[0] != 0x0){free(btp.NT[0]);}
 if(btp.NT[1] != NULL && btp.NT[1] != 0x0){free(btp.NT[1]);}
 //Print a warning if the void pointer doesn't seem to be empty
 if(btp.VP!= NULL && btp.VP != 0x0){fprintf(stderr,"WARNING: deallocating bond_type structure with non-NULL void pointer!\n");}
 return ;
}
void deallocateAngleType(angle_type atp){
 //Free up the description
 if(atp.desc != NULL && atp.desc != 0x0){free(atp.desc);}
 //Free up the first and second atom types in bond
 if(atp.NT[0] != NULL && atp.NT[0] != 0x0){free(atp.NT[0]);}
 if(atp.NT[1] != NULL && atp.NT[1] != 0x0){free(atp.NT[1]);}
 if(atp.NT[2] != NULL && atp.NT[2] != 0x0){free(atp.NT[2]);}
 //Print a warning if the void pointer doesn't seem to be empty
 if(atp.VP!= NULL && atp.VP != 0x0){fprintf(stderr,"WARNING: deallocating angle_type structure with non-NULL void pointer!\n");}
 return ;
void deallocateTorsionType(torsion_type ttp){
 //Free up the description
 if(ttp.desc != NULL && ttp.desc != 0x0){free(ttp.desc);}
 //Free up the first and second atom types in bond
 if(ttp.NT[0] != NULL && ttp.NT[0] != 0x0){free(ttp.NT[0]);}
 if(ttp.NT[1] != NULL && ttp.NT[1] != 0x0){free(ttp.NT[1]);}
 if(ttp.NT[2] != NULL && ttp.NT[2] != 0x0){free(ttp.NT[2]);}
 if(ttp.NT[3] != NULL && ttp.NT[3] != 0x0){free(ttp.NT[3]);}
 //Free the n terms
 if(ttp.k != NULL && ttp.k != 0x0){free(ttp.k);} 
 if(ttp.N != NULL && ttp.N != 0x0){free(ttp.N);} 
 if(ttp.P != NULL && ttp.P != 0x0){free(ttp.P);} 
 //Print a warning if the void pointer doesn't seem to be empty
 if(ttp.VP!= NULL && ttp.VP != 0x0){fprintf(stderr,"WARNING: deallocating torsion_type structure with non-NULL void pointer!\n");}
 return ;
}

void deallocateAtype(atype atp){
 int i;
 if(atp.N != NULL && atp.N != 0x0){free(atp.N);}
 if(atp.NT != NULL && atp.NT != 0x0){free(atp.NT);}
 if(atp.desc != NULL && atp.desc != 0x0){free(atp.desc);}
 if(atp.bo != NULL && atp.bo != 0x0){free(atp.bo);}
 if(atp.ch != NULL && atp.ch != 0x0){free(atp.ch);}
 if(atp.R != NULL && atp.R != 0x0){free(atp.R);}
 if(atp.SC != NULL && atp.SC != 0x0){free(atp.SC);}
// free up the convenience indices & pointers
 if(atp.iBT != NULL && atp.iBT != 0x0){free(atp.iBT);}
 if(atp.BT != NULL && atp.BT != 0x0){free(atp.BT);}
 if(atp.iHBT != NULL && atp.iHBT != 0x0){free(atp.iHBT);}
 if(atp.HBT != NULL && atp.HBT != 0x0){free(atp.HBT);}
 if(atp.iNBT != NULL && atp.iNBT != 0x0){free(atp.iNBT);}
 if(atp.NBT != NULL && atp.NBT != 0x0){free(atp.NBT);}
 if(atp.iANT != NULL && atp.iANT != 0x0){free(atp.iANT);}
 if(atp.ANT != NULL && atp.ANT != 0x0){free(atp.ANT);}
 if(atp.iHANT != NULL && atp.iHANT != 0x0){free(atp.iHANT);}
 if(atp.HANT != NULL && atp.HANT != 0x0){free(atp.HANT);}
 if(atp.iNANT != NULL && atp.iNANT != 0x0){free(atp.iNANT);}
 if(atp.NANT != NULL && atp.NANT != 0x0){free(atp.NANT);}
 if(atp.iTRT != NULL && atp.iTRT != 0x0){free(atp.iTRT);}
 if(atp.TRT != NULL && atp.TRT != 0x0){free(atp.TRT);}
 if(atp.iHTRT != NULL && atp.iHTRT != 0x0){free(atp.iHTRT);}
 if(atp.HTRT != NULL && atp.HTRT != 0x0){free(atp.HTRT);}
 if(atp.iNTRT != NULL && atp.iNTRT != 0x0){free(atp.iNTRT);}
 if(atp.NTRT != NULL && atp.NTRT != 0x0){free(atp.NTRT);}
 //Free n terms
 for(i=0;i<atp.nR;i++){ if(atp.RD[i] != NULL && atp.RD[i] != 0x0){free(atp.RD[i]);} } 
 //Print a warning if the void pointer doesn't seem to be empty
 if(atp.VP!= NULL && atp.VP != 0x0){fprintf(stderr,"WARNING: deallocating atom_type structure with non-NULL void pointer!\n");}
 return ;
} 
void deallocateRtype(rtype rtp){
 if(rtp.ac != NULL && rtp.ac != 0x0){free(rtp.ac);} 
 //Print a warning if the void pointer doesn't seem to be empty
 if(rtp.VP!= NULL && rtp.VP != 0x0){fprintf(stderr,"WARNING: deallocating rtype structure with non-NULL void pointer!\n");}
 return ;
}
void deallocateMtype(mtype mtp){
 //Print a warning if the void pointer doesn't seem to be empty
 if(mtp.VP!= NULL && mtp.VP != 0x0){fprintf(stderr,"WARNING: deallocating mtype structure with non-NULL void pointer!\n");}
 return ;
} 
void deallocateTypes(types tps){
 int i;
 //deallocate substructures
	for(i=0;i<tps.na;i++) deallocateAtype(tps.a[i]);
	for(i=0;i<tps.nr;i++) deallocateRtype(tps.r[i]);
	for(i=0;i<tps.nm;i++) deallocateMtype(tps.m[i]);
 //Print a warning if the void pointer doesn't seem to be empty
 if(foo.VP!= NULL && foo.VP != 0x0){fprintf(stderr,"WARNING: deallocating foo structure with non-NULL void pointer!\n");}
 return ;
}

void deallocateParameterSet(parameter_set ps){
 int i;
 //deallocate substructures
	for(i=0;i<tps.nAT;i++) deallocateAtype(tps.AT[i]);
	for(i=0;i<tps.nRT;i++) deallocateRtype(tps.RT[i]);
	for(i=0;i<tps.nMT;i++) deallocateMtype(tps.MT[i]); 
	for(i=0;i<tps.nBT;i++) deallocateBondType(tps.BT[i]);
	for(i=0;i<tps.nHBT;i++) deallocateBondType(tps.HBT[i]);
	for(i=0;i<tps.nNBT;i++) deallocateBondType(tps.NBT[i]);
	for(i=0;i<tps.nANT;i++) deallocateAngleType(tps.ANT[i]);
	for(i=0;i<tps.nTRT;i++) deallocateTorsionType(tps.TRT[i]);
 return ;
}

/********** structures from molecules.h ****************/

void deallocateRingEnsindex(ring_ensindex re){
 if(re.P != NULL && re.P != 0x0){free(re.P);} 
 return ;
}
typedef struct {
	int nP; ///< number of positions in the ring
	ensindex *P; ///< the nP relevant positions
	int nin,*in; ///< reference integers in *P for other ring members with incoming bonds
	int nout,*out; ///< reference integers in *P for other ring members with outgoing bonds
} ring_ensindex; ///< Structure holding ensemble indices for a ring, plus maybe other info
void deallocatefoo(foo bar){
 int i;
	if( != NULL &&  != 0x0){free();} 
	for(i=0;i< ;i++){
		if( != NULL &&  != 0x0){free();} 
		}
 //Print a warning if the void pointer doesn't seem to be empty
 if(foo.VP!= NULL && foo.VP != 0x0){fprintf(stderr,"WARNING: deallocating foo structure with non-NULL void pointer!\n");}
 return ;
}
typedef struct { 
	int na, *ai; ///< na atom indices and the na indices 
} residue_tree_index; ///< for storing indices within a structure like E[ei].m[mi].r[ri].a[ai].i
void deallocatefoo(foo bar){
 int i;
	if( != NULL &&  != 0x0){free();} 
	for(i=0;i< ;i++){
		if( != NULL &&  != 0x0){free();} 
		}
 //Print a warning if the void pointer doesn't seem to be empty
 if(foo.VP!= NULL && foo.VP != 0x0){fprintf(stderr,"WARNING: deallocating foo structure with non-NULL void pointer!\n");}
 return ;
}
typedef struct { 
	int nr, *ri; ///< nr residue indices and the nr indices
	residue_tree_index *r; ///< nr residue_tree_index structures
} molecule_tree_index; ///< for storing indices within a structure like E[ei].m[mi].r[ri].a[ai].i
void deallocatefoo(foo bar){
 int i;
	if( != NULL &&  != 0x0){free();} 
	for(i=0;i< ;i++){
		if( != NULL &&  != 0x0){free();} 
		}
 //Print a warning if the void pointer doesn't seem to be empty
 if(foo.VP!= NULL && foo.VP != 0x0){fprintf(stderr,"WARNING: deallocating foo structure with non-NULL void pointer!\n");}
 return ;
}
typedef struct { 
	int nm, *mi; ///< nm molecule indices and the nm indices
	molecule_tree_index *m; ///< nm molecule_tree_index structures
} ensemble_tree_index; ///< for storing indices within a structure like E[ei].m[mi].r[ri].a[ai].i

void deallocatefoo(foo bar){
 int i;
	if( != NULL &&  != 0x0){free();} 
	for(i=0;i< ;i++){
		if( != NULL &&  != 0x0){free();} 
		}
 //Print a warning if the void pointer doesn't seem to be empty
 if(foo.VP!= NULL && foo.VP != 0x0){fprintf(stderr,"WARNING: deallocating foo structure with non-NULL void pointer!\n");}
 return ;
}
typedef struct {
	int posneg; ///< Is this a positive (1), negative (-1) or brand new (0) selection set?
	int nmN,nmi,nmn,*mn,*mi; ///< # of molecule names, numbers & indicies, nmn numbers and nmi indices
	char **mN; ///< nmN names
	int nrN,nri,nrn,*rn,*ri; ///< # of residue names, numbers & indicies, nrn numbers and nri indices
	char **rN; ///< nrN names
	int naN,nai,nan,*an,*ai; ///< # of atom names, numbers & indicies, nan numbers and nai indices
	char **aN; ///< naN names
} moiety_selection; 

void deallocatefoo(foo bar){
 int i;
	if( != NULL &&  != 0x0){free();} 
	for(i=0;i< ;i++){
		if( != NULL &&  != 0x0){free();} 
		}
 //Print a warning if the void pointer doesn't seem to be empty
 if(foo.VP!= NULL && foo.VP != 0x0){fprintf(stderr,"WARNING: deallocating foo structure with non-NULL void pointer!\n");}
 return ;
}
typedef struct { 
	ensindex s; ///< "source" -- index to first atom in bond
	ensindex t; ///< "target" -- index to the other atom in the bond
	double o; ///< order of bond
	bond_type *typ; ///< pointer to type of bond
	int i; ///< index -- for example as alternative to *typ
	char *D; ///< free-form description
} bond; ///< a generic bond structure for connections of any type
typedef struct {
	int n; // number of bonds
	bond *b; // n of these
} bondset; // set of consecutive bonds


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
