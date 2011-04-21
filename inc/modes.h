/** \file  modes.h
\brief  Header file for vibrational modes.

Begun 20070519 by Lachele Foley
This file contains structures to be used for descriptions of
vibrational modes.  Most of the definitions should be obvious
to most folks who are concerned with vibrations.  But, a few,
like the "extent of" descriptors, will not be so evident.  
Eventually, I'll document how they are determined -- they
basically boil down to normalized dot products, but the exact
definitions might be more subtle. */

#if !defined(GLYLIB_MODES)
#define GLYLIB_MODES
//#include "molecules.h"

/** \addtogroup VIBRATIONS
 * @{
 */
typedef struct {
	int i; ///< array index
	int j; ///< internal index 
	int k; ///< some other index
} intset; ///< intended to use in assignment arrays for indexing to mode details

typedef struct {
	int ct; ///< mode class type (see docs)
// Class types are indexed by number, e.g.:
//	0	stretch  
//	1	bend 
//	2	torsion 
//	3	ringmotion
	int ce; ///< class entry ... 
	int i ; ///< other index
// the index for of a certain type (e.g., stretch) in the vibmode array
	int C; ///< complexity for this mode entry
	double w; ///< weight (intensity, extent, etc.) for this mode
} vibindex;

typedef struct{
	int i; ///< general index
	int nv; ///< number of va's
	vibindex *vi; ///< indices for mode details
	molindex mi; ///< info about the atom
} vibaddr; ///< for holding addresses in a vibmode array

typedef struct {
	int i; ///< index
	char *Desc; ///< free-form description for mode 
	double I; ///< intensity of the motion (user-defined...)
	double e; ///< extent of stretch in motion
	molindex  A,B;
	coord_3D c; ///< e.g., for describing motions in i,j,k 
} stretch; ///< description of stretch

typedef struct {
	int i; ///< index
	char *Desc; ///< free-form description for mode 
	double I; ///< intensity of the motion (user-defined...)
	double e; ///< extent of angle bend in motion
	molindex  A,B,C; 
	coord_3D c; ///< e.g., for describing motions in i,j,k 
} bend; 

typedef struct {
	int i; ///< index
	char *Desc; ///< free-form description for mode 
	double I; ///< intensity of the motion (user-defined...)
	double e; ///< extent of torsion in motion
	molindex  A,B,C;
	coord_3D c; ///< e.g., for describing motions in i,j,k 
} torsion; ///< description of torsion

typedef struct {
	int i; ///< index
	char *Desc; ///< free-form description for mode
	char *RingDesc; ///< Description of ring type, e.g. phenyl, 0GB, etc.
	double I; ///< intensity of the motion (user-defined...)
	double e; ///< extent of ring-portion of overall motion
	double pvp; ///< motion parallel (~0) or perpendicular (~1) to plane
	double par; ///< breathing-type or side-to-side motion (if parallel)
	double sym; ///< symmetric (~1) or antisymmetric (~0) (par or perp)
	int symZ,symXY; ///< symmetries in the Z direction and XY plane
	int na; ///< number of main atoms involved
	molindex  *mi; ///<na of these
	int nc; ///< number of vectormag_3D sets 
	vectormag_3D *c; ///< nc of these (e.g., direction intensities of atoms)
} ringmotion; ///< description of ring-related motion

typedef struct {
	int i,j,k; ///< indices -- molecule number, file entry number, etc.
	double f; ///< frequency for this mode (usually in wavenumbers)
	double I; ///< intensity of mode, user-defined
	double IR; ///< IR intensity of the mode (absorption, preferrably)
	double RA; ///< RAMAN intensity of the mode (absorption, preferrably)
	char *Desc; ///< free-form description
	int C; ///< complexity factor (number of different motions represented)
	int ns; ///< number of stretches in mode
	int nb; ///< number of angle bends in mode
	int nt; ///< number of torsions in mode
	int nr; ///< number of ring-related motions in mode
	stretch *s; ///< stretches (ns of these)
	bend *b; ///< number of angle bends (na of these)
	torsion *t; ///< torsions (nt of these)
	ringmotion *r; ///< ring-related motions (nr of these)
	double te; ///< extent of translational character in mode
	double re; ///< extent of rotational character in mode
	char ismoltors; ///< 'y' if the whole molecule is torsioning, 'n' if not
	vectormag_3D tv, rv, mvs; ///< translation vec, rotation vec, sum of motion vectors
	double maxmag; ///< maximum vector magnitude in native set
	double modemax; ///< maximum s-a-t score in mode set (used to scale vectors to max of 1)
	double maxmodemag; ///< vector magnitude for the modemax vector
	double maxKEmag; ///< for the larges s-a-t assignment, the KE along that assignment
	double TOTmag; ///< Total magnitude of all motions, all molecules
	double TRmag; ///< Total translation-removed magnitude, all molecules
	double MRmag; ///< Total molecule-rotation-removed magnitude, all molecules
	double RRmag; ///< Total residue-rotation-removed magnitude, all molecules
} vibmode; ///< description of a vibrational mode 


typedef struct { // human-readable-assignment structure for two equivalent atoms
	int C; ///<complexity
	int ns2,nb2,nt2; ///< for atom list with unknown symmetry
	intset *s2,*b2,*t2; 
	int ns2s;
	intset *s2s; ///< stretches with two equivalent atoms -- symmetric
	int ns2a;
	intset *s2a; ///< stretches with two equivalent atoms -- asymmetric
	int nb2s;
	intset *b2s; ///< stretches with two equivalent atoms -- symmetric
	int nb2a;
	intset *b2a; ///< stretches with two equivalent atoms -- asymmetric 
	int nt2s;
	intset *t2s; ///< stretches with two equivalent atoms -- symmetric
	int nt2a;
	intset *t2a; ///< stretches with two equivalent atoms -- asymmetric
	double s2se,s2ae,b2se,b2ae,t2se,t2ae; ///< these are "extents" for each of the above motions
} twoassign; ///< description of a vibrational mode 
typedef struct { // human-readable-assignment structure for two ring-related motions
	int C; ///<complexity
	int nopsA;
	intset *opsA; ///< Out of plane stretch (axial atoms), single atom info
//	-- these point to ringmotion structures relating to the atoms above
	int nopss;
	intset  *opss; ///< Out of plane stretch (axial atoms) -- symmetric
	int nopsa;
	intset  *opsa; ///< Out of plane stretch (axial atoms) -- asymmetric
	double opsSe,opsPe,opsDe; ///< these are "extents" for each of the above motions
	int nopdA;
	intset  *opdA; ///< Out of plane distortion (ring atoms), single atom info
//	-- these point to ringmotion structures relating to the atoms above
	int nopdA1G;
	intset  *opdA1G; ///< Out of plane distortion (ring atoms) -- A-1-G symmetry
	int nopdOth;
	intset  *opdOth; ///< Out of plane distortion (ring atoms) -- some other symmetry
	double opdse,opdae; ///< these are "extents" for each of the above motions
	int nrbA;
	intset  *rbA; ///< ring breathing motion, single atom info
//	-- these point to ringmotion structures relating to the atoms above
	int nrbs;
	intset  *rbs; ///< ring breathing motion, symmetric
	int nrba;
	intset  *rba; ///< ring breathing motion, asymmetric
// REGARDING O and P assignments:
// 	Atom A is exocyclic and bound to endocyclic B
//	R is the vector for the radius from ring center to B
//	T is the tangent at B, always oriented the same, e.g. (counter)clockwise
//	M is the cross product of R and T (M points perpendicular to the ring plane)
// If the motion vector is orthogonal to M, it is type O (moves along ring perimeter)
// If the motion vector is parallel to M, it is type P (moves _|_ to ring perimeter)
	int nroA;
	intset  *roA; ///< ring-related O motion (not stretch), single atom info
//	-- these point to ringmotion structures relating to the atoms above
	int nros;
	intset  *ros; ///< ring-related O motion, symmetric
	int nroa;
	intset  *roa; ///< ring-related O motion, asymmetric
	double rbAe,rbse,rbae,roAe,rose,roae;  ///< these are "extents" for each of the above motions
	int nrpA;
	intset  *rpA; ///< ring-related P motion, (not breathing) single atom info
//	-- these point to ringmotion structures relating to the atoms above
	int nrps;
	intset  *rps; ///< ring-related P motion, symmetric
	int nrpa;
	intset  *rpa; ///< ring-related P motion, asymmetric
	int nrso;
	intset  *rso; ///< other ring-related stretches
	double rpAe,rpse,rpae,rsoe;  ///< these are "extents" for each of the above motions
} ringassign; ///< description of a vibrational mode 
typedef struct { // human-readable-assignment structure
	int i; ///< index
	int n; ///< index (intended as atom index for per-atom info)
	int loc; ///< location index (endocyclic, exocyclic H, etc.)
	char *Desc; ///< free-form description
	int C; ///< complexity = number of motion types present 
// the following pairs are the number per each and their index in an external vibmode array
//	Simple motions
	int ns;
	intset *s; ///< plain stretches in mode 
	int nb;
	intset *b; ///< plain angle bends in mode
	int nt;
	intset *t; ///< plain torsions in mode
	int nw;
	intset *w; ///< number of plain wags, pointers to locations
	int nrk;
	intset *rk; ///< number of plain rocks, pointers to locations
	double se,be,te,we,rke; ///< these are "extents" for each of the above motions
// 	Motions of two coupled atoms and motions involving rings
	int nr, ntwo; ///< numbers of rings or pairs of coupled atoms
	ringassign *r; ///< ring assignment structures
	twoassign *two; ///< structures for coupled motions
} assignment;

typedef struct{
	int i,j,k; ///< indices
	molindex mi; ///< where this atom is
	double a,b,c; ///< numbers of interest
	char *Desc; ///< description 
	//intset aset,bset; // oh, why not...
} assn_brief;

typedef struct { // info tied to atom number
	int i,j; ///< index (to an assignment array, perhaps)
	double I; ///< normalized motion intensity for this atom (regular intensity is elsewhere)
	int ns, *s; ///< number of significant stretches, indices for those in VS
	int nb, *b; ///< number of significant bends, indices for those in VS
	int nt, *t; ///< number of significant torsions, indices for those in VS
	int loc; // the location of this atom
} atommode;

void get_transition_intensity(molecule *M, int xl, int vl, vibmode *v, atype *ATYPE, FILE *F);
/** @}*/
#endif
