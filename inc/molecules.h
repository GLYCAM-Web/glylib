/** \file molecules.h
	Header file for molecule structures 

File begun May 2007 by Lachele Foley and modified continually ever since
	by Lachele Foley and others, including Mike Tessier. */

#if !defined(GLYLIB_MOLECULES)
#define GLYLIB_MOLECULES
#include <geometries.h>
#include <parameter_sets.h>

char *ATYPESFILE,*RTYPESFILE,*MTYPESFILE; ///< locations of type databases,
	// to be set by the program, but globally visible

// these now in geometries.h
// coord_3D; // cartesian, polar, direction cosines, etc.
// vectormag_3D; // same as "plane" but with easier naming
// plane; // standard plane, Ax+By+Cz+D=0

// these now in parameter_sets.h, and greatly expanded
// atype; // atom types, with other info
// rtype; // residue types, with other info
// mtype; // molecule types, with oter info
// types; // superstructure for typing information

typedef struct {
	int i; ///< general index
	int m; ///< molecule index
	int r; ///< residue index
	int a; ///< atom index
} molindex; ///< Index to describe position in a molecule
typedef struct {
	int i; ///< general index
	int E; ///< ensemble index
	int A; ///< assembly index
	int m; ///< molecule index
	int r; ///< residue index
	int a; ///< atom index
} ensindex; ///< Index to describe position in an ensemble
typedef struct {
	int nP; ///< number of positions in the ring
	ensindex *P; ///< the nP relevant positions
	int nin,*in; ///< reference integers in *P for other ring members with incoming bonds
	int nout,*out; ///< reference integers in *P for other ring members with outgoing bonds
} ring_ensindex; ///< Structure holding ensemble indices for a ring, plus maybe other info

typedef struct { 
	int na, *ai; ///< na atom indices and the na indices 
} residue_tree_index; ///< for storing indices within a structure like E[ei].m[mi].r[ri].a[ai].i
typedef struct { 
	int nr, *ri; ///< nr residue indices and the nr indices
	residue_tree_index *r; ///< nr residue_tree_index structures
} molecule_tree_index; ///< for storing indices within a structure like E[ei].m[mi].r[ri].a[ai].i
typedef struct { 
	int nm, *mi; ///< nm molecule indices and the nm indices
	molecule_tree_index *m; ///< nm molecule_tree_index structures
} ensemble_tree_index; ///< for storing indices within a structure like E[ei].m[mi].r[ri].a[ai].i

typedef struct {
	int posneg; ///< Is this a positive (1), negative (-1) or brand new (0) selection set?
	int nmN,nmi,nmn,*mn,*mi; ///< # of molecule names, numbers & indicies, nmn numbers and nmi indices
	char **mN; ///< nmN names
	int nrN,nri,nrn,*rn,*ri; ///< # of residue names, numbers & indicies, nrn numbers and nri indices
	char **rN; ///< nrN names
	int naN,nai,nan,*an,*ai; ///< # of atom names, numbers & indicies, nan numbers and nai indices
	char **aN; ///< naN names
} moiety_selection; 

/* 20100127 BLFoley: I'm changing the bond structure.  With luck, this will 
only break the vibrations program I'm rewriting... */
/// for simple, all-one-residue bonding (for early programs)
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

/* The following is a set of functions intended primarily to facilitate 
debugging.  But, they could easily be used for verbose data output within
the program.  They all follow the same basic format:

print_X(X *x, int level)

where 	X = one of the structures defined above (e.g., molecule, residue, 
		bondset, atom, bond, atype, plane, vectormag_3D, coord_3D).
	*x = pointer to your structure
	level = the depth to print.  This is only available on structures 
		likely to contain arrays of other structures.  level=0 only
		prints the top-level information.  level=1 prints one level 
		down of other structures, lists, etc. level=2 prints two 
		levels down, etc.   The value of level can be larger than 
		the depth available -- the function will stop when it runs
		out of structures to print.  
*/
int NUMAT; 
void print_molecule(molecule*,int),print_residue(residue*,int);
void print_bondset(bondset*,int),print_atom(atom*,int),print_bond(bond*);
void print_atype(atype*,int),print_plane(plane*);
void print_vectormag_3D(vectormag_3D*),print_coord_3D(coord_3D*);
//
void dprint_molecule(molecule*,int),dprint_residue(residue*,int);
void dprint_bondset(bondset*,int),dprint_atom(atom*,int),dprint_bond(bond*);
void dprint_atype(atype*,int),dprint_plane(plane*);
void dprint_vectormag_3D(vectormag_3D*),dprint_coord_3D(coord_3D*);
//
void dXprint_molecule(molecule*,int),dXprint_residue(residue*,int);
void dXprint_bondset(bondset*,int),dXprint_atom(atom*,int),dXprint_bond(bond*);
void dXprint_atype(atype*,int),dXprint_plane(plane*);
void dXprint_vectormag_3D(vectormag_3D*),dXprint_coord_3D(coord_3D*);

/* the following are geometric operations for the structures defined above 
	where int is xl ("x" (coordinate) location):
		xl refers to position in the atom structure
			xl=-1 means use the atom.x structure
			xl=0,1,2,etc means use the atom.xa[xl] or
				the atom.v[xl] structure
	where int,int is "xl,vl" it refers to coord and vector location
		(in that order).
	where int is n
		n indicates the number of blocks in the preceding pointer

   All of the following that involve vectors will set the magnitude of 
	the resulting vector without relying on the entry in the input 
	vector.  So, for example, there is no need to get_magnitude 
	before normalizing or multiplying by a scalar.
*/
/* this moves all the coords in molecule so that the center of mass
for the molecule is at the origin (right now, though, it only works
if the molecule has only one residue...*/
void translate_to_COM(molecule *,atype *,int); // int is xl
	// this assumes the initial coords are in atom.x
//void translate_zero_to_COM(molecule *,atype *,int,int); 
	// int#1 is xs, the location of the source coords
	// int#2 is xt, the location of the translated coords
void translate_zero_to_coord_M(molecule *,int,int,coord_3D); 
	// int#1 is xs, the location of the source coords
	// int#2 is xt, the location of the translated coords
void assign_residue_COM(residue *r,atype *ATYPE);
void assign_molecule_COM(molecule *m,atype *ATYPE);
void assign_assembly_COM(assembly *a,atype *ATYPE);
void assign_ensemble_COM(ensemble *e,atype *ATYPE);
coord_3D get_residue_COM(residue *r,atype *ATYPE,int xs);
coord_3D get_molecule_COM(molecule *m,atype *ATYPE,int xs);
coord_3D get_assembly_COM(assembly *a,atype *ATYPE,int xs);
coord_3D get_ensemble_COM(ensemble *e,atype *ATYPE,int xs);

/* the following functions rotate the coordinate list such that the
vector given by vectormag_3D points along the Z (Y, X, other) axis.
The orientation of the orthogonal axes is arbitrary. */
void rotate_vector_to_Z_M(molecule*,int,int,int,int,vectormag_3D); 
	// The four integers are:
	// int#1 location of source coordinates (-1=x, 0,1,2,etc=xa[xs]
	// int#2 location of rotated coordinates (-1=x, 0,1,2,etc=xa[xs]
	// int#3 location of source vectors (-1=do not calc, 0,1,2,etc=v[xs]
	// int#4 location of rotated vectors (-1=do not calc, 0,1,2,etc=v[xs]
void rotate_vector_to_Z_list(coord_3D*,int,vectormag_3D); // int is n
//void rotate_vector_to_Y_M(molecule*,int,vectormag_3D); // int is xl
//void rotate_vector_to_Y_list(coord_3D*,int,vectormag_3D); // int is n
//void rotate_vector_to_X_M(molecule*,int,vectormag_3D); // int is xl
//void rotate_vector_to_X_list(coord_3D*,int,vectormag_3D); // int is n
/* in these, the firt vector is rotated onto the second vector */
//void rotate_vector_to_V_M(molecule*,int,vectormag_3D,vectormag_3D); // int is xl
//void rotate_vector_to_V_list(coord_3D*,int,vectormag_3D,vectormag_3D); // int is n

coord_3D get_geometric_center(coord_3D *,int); // int is n
plane get_plane(coord_3D,coord_3D,coord_3D);
vectormag_3D normalize_vec(vectormag_3D);
void normalize_molecule_vectors(molecule *,int,int);
	// int #1 -- location of source vector
	// int #2 -- where to write vectors normalized to one (made unit vectors) 
vectormag_3D scalarmult_vec(vectormag_3D,double);
vectormag_3D add_vec(vectormag_3D,vectormag_3D); // add two vectors
vectormag_3D subtract_vec(vectormag_3D,vectormag_3D); // subtract second vector from first
coord_3D scalarmult_coord(coord_3D,double); // add two vectors
coord_3D add_coord(coord_3D,coord_3D); // add two vectors
coord_3D subtract_coord(coord_3D,coord_3D); // subtract second vector from first
vectormag_3D get_crossprod(vectormag_3D, vectormag_3D); // returns the cross product
double get_dotprod(vectormag_3D, vectormag_3D); // returns dot product of two vectors
/* the following essentially returns the cosine of the angle between two vectors */
//double get_dotprodN(vectormag_3D, vectormag_3D); // returns dot prod, but normalizes vecs first
double get_magnitude(vectormag_3D); // calculates vector magnitude for "d" in structure
void shift_molecule_atoms_by_vector_scale(molecule *m,int xs,int xt,int vr, double scale);
vectormag_3D get_molecule_point_charge_dipole(molecule *m,int xs,int chgsrc, atype *AT);

vectormag_3D zero_vec(); // zeros a vector
coord_3D zero_coord(); // zeros a coordinate set
coord_3D vec_to_coord(vectormag_3D); // turns a vector into a coordinate set
vectormag_3D coord_to_vec(coord_3D); // turns a coordinate set into a vector
void rollMolecule(molecule*,double); //Rotates about x-axis using radians
void pitchMolecule(molecule*,double);//Rotates about y-axis using radians
void yawMolecule(molecule*,double);  //Rotates about z-axis using radians
void rollAssembly(assembly*,double); //Rotates about x-axis using radians
void pitchAssembly(assembly*,double);//Rotates about y-axis using radians
void yawAssembly(assembly*,double);  //Rotates about z-axis using radians
void outputMolPDB(molecule*,char*);  //Writes a pdb using a given molecule
void outputAsmblPDB(assembly*,char*);//Writes a pdb using a given assembly
char* spacing(int,int); //Used by output functions to determine spacing

// RMS between coordinate sets xs and xt
double get_alt_rms_res(residue *r, int xs, int xt); // per residue
double get_alt_rms_mol(molecule *m, int xs, int xt); // per molecule

atype *ATYPE_init(); // initialize atom types -- for a few old programs only
// use void load_atypes(fileset FT, types *T); instead (declarations.h)
//dockinfo *load_dlg_mol(fileset F,atype *AT); also in declarations.h

void initialize_coord_3D(coord_3D *c);
void initialize_vectormag_3D(vectormag_3D *v);
void initialize_plane(plane *p);
void initialize_atype(atype *at);
void initialize_rtype(rtype *rt);
void initialize_mtype(mtype *mt);
void initialize_types(types *t);
void initialize_molindex(molindex *mi);
void initialize_ensindex(ensindex *ei);
void initialize_bond(bond *b);
void initialize_bondset(bondset *bs);
void initialize_molbond(molbond *mb);
void initialize_molbondset(molbondset *mbs);
void initialize_atom(atom *a);
void initialize_residue(residue *r);
void initialize_molecule(molecule *m);
void initialize_dockinfo(dockinfo *di);
void initialize_assembly(assembly *A);
void initialize_ensemble(ensemble *E);

//Functions that recursively free memory in themselves, with the exception of parameters
void deallocateBondType(bond_type btp);
void deallocateMolbond(molbond mlb);
void deallocateConnectionTree(connection_tree con);
void deallocateBond(bond bnd);
void deallocateBondset(bondset bst);
void deallocateAtom(atom atm);
void deallocateResidue(residue res);
void deallocateMolecule(molecule mol);

//Functions that add or remove structures from other structures
void add_assembly_to_ensemble(
        assembly *A, ///< pointer to the assembly being added (SEE DOCS)
        ensemble *E ///< pointer to the ensemble being grown
        );

#endif
