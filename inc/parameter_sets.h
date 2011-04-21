/** \file parameter_sets.h 
\brief Structures for adding parameters to ensembles and assemblies.

Begun on 20080606 by BLFoley

This file is loaded with the molecules.h header file by default.
	It should only be loaded explicitly in a program if, for some reason,
	you are using it separately. 

The parameter_set structure will eventually contain pointers for
just about any type of data typically used in simulations.  As usual,
there is also a void pointer for situations where information is needed
that is not already included in the structure.  In any case, this is
likely to be one whopping big structure, hence the separate file. */

#if !defined(GLYLIB_PARAMETERS)
#define GLYLIB_PARAMETERS

/** \addtogroup PARAMETERS
 * @{
 */
/************ structures for descriptions of molecules ***********************/
typedef struct {
	char RS,DL,pm; ///< for R/S & D/L type chirality & p/m type optical activity
	char UD, AB; ///< for up/down and alpha/beta descriptors
	int niso; ///< number of other descriptors
	char **iso; ///< the niso other descriptions of isomerism
	char *ELGEOM,*SPGEOM; ///< electronic (sp2, sp3, etc.) and spatial (tetrahedral, etc.) geometry
} chirality_description;


/************ sub structures for the parameter_set structure *****************/
typedef struct{ // bond_type
	char **NT; // first and second atom types in bond (match NT in atype)
	char *desc; // description 
	double o; // bond order for this type
	double k; // force constant 
	double l; // equilibrium length
	double kx; // anharmonicity 
	double LJD_612,LJZ_612,LJD_1012,LJZ_1012; // Lennard-Jones params (D=depth Z=zero)
	double LJ6_612,LJ12_612,LJ10_1012,LJ12_1012; // Lennard-Jones params (alternate)
	int nVP ; // number void pointers
	void *VP; // the void pointers
} bond_type; // bond types -- for describing bonds
typedef struct { // angle_type
	char **NT; // the three atom types (match NT in atype)
	char *desc; // description 
	double k; // force constant (angular)
	double l; // equilibrium angle
	double kx; // anharmonicity 
	int nVP ; // number void pointers
	void *VP; // the void pointers
} angle_type; // for describing angle parameters
typedef struct { // torsion_type
	char **NT; // the four atom types (match NT in atype)
	char *desc; // description 
	int n; // number of terms in this dependence
	double *k; // force constant (angular) (n of these...)
	double *N; // periodicity (n of these...)
	double *P; // phase (n of these...)
	int nVP ; // number void pointers
	void *VP; // the void pointers 
} torsion_type; // for describing torsion parameters
typedef struct { // atype -- for atom types
	int n; // atomic number
	char *N; // element name
	char *NT; // name for element type
	double m; // mass of this atom type 
	char *desc; // brief free-form description field
	int nb; // number of typical bonds for this atom type
	double *bo; // typical bond orders (nb of these)
	int nlp; // number of typical lone pairs (to check geometry sanity)
	double LJD_612,LJZ_612,LJD_1012,LJZ_1012; // Lennard-Jones params (D=depth Z=zero)
	int nch; // number of charges to associate with this atom type
	int *ch; // the charges
	int nR; // number of radii 
	char **RD; // descriptions of the nR radii
	double *R; // radii (e.g., van der Waals, Poisson-Boltzmann, etc.)
	int nSC; 
	double *SC; // some other constant (e.g., screening) associated with R
// Pointers for convenience
	int nBT; // number of bond types for this atom
	int *iBT; // pointers into a bond_type array where likely info is stored (to reduce search time) 
	bond_type **BT; // pointers into a bond_type array where likely info is stored (to reduce search time) 
	int nHBT; // number of H-bond types for this atom
	int *iHBT; // pointers into a H-bond_type array where likely info is stored (to reduce search time) 
	bond_type **HBT; // pointers into a H-bond_type array where likely info is stored (to reduce search time) 
	int nNBT; // number of non-bond types for this atom
	int *iNBT; // pointers into a non-bond_type array where likely info is stored (to reduce search time) 
	bond_type **NBT; // pointers into a non-bond_type array where likely info is stored (to reduce search time) 
	int nANT; // number of angle types for this atom
	int *iANT; // pointers into a angle_type array where likely info is stored (to reduce search time) 
	angle_type **ANT; // pointers into a angle_type array where likely info is stored (to reduce search time) 
	int nHANT; // number of H-bond angle types for this atom
	int *iHANT; // pointers into a H-bond angle_type array where likely info is stored (to reduce search time) 
	angle_type **HANT; // pointers into a H-bond angle_type array where likely info is stored (to reduce search time) 
	int nNANT; // number of non-bond angle_types for this atom
	int *iNANT; // pointers into a non-bond angle_type array where likely info is stored (to reduce search time) 
	angle_type **NANT; // pointers into a non-bond angle_type array where likely info is stored (to reduce search time) 
	int nTRT; // number of torsion types for this atom
	int *iTRT; // pointers into a torsion_type array where likely info is stored (to reduce search time) 
	torsion_type **TRT; // pointers into a torsion_type array where likely info is stored (to reduce search time) 
	int nHTRT; // number of H-bond torsion types for this atom
	int *iHTRT; // pointers into a H-bond_torsion_type array where likely info is stored (to reduce search time) 
	torsion_type **HTRT; // pointers into a H-bond_torsion_type array where likely info is stored (to reduce search time) 
	int nNTRT; // number of non-bond torsion_types for this atom
	int *iNTRT; // pointers into a non-bond torsion_type array where likely info is stored (to reduce search time) 
	torsion_type **NTRT; // pointers into a non-bond torsion_type array where likely info is stored (to reduce search time) 
// The void pointers for expansion
	int nVP ; // number void pointers
	void *VP; // the void pointers
} atype; // atom types, with other info 
typedef struct { // rtype -- for residue types
	int c; // main class (amino acid, glycan, solvent, etc.)
	int nac, *ac; // number alternate classes, those classes
	int nVP; // number of other information
	void *VP; // pointers to the nVP holders of other information
} rtype; // residue types, with other info 
typedef struct { // mtype -- for molecule types
	int c; // type class (amino acid, glycan, solvent, etc.)
	int nVP; // number of other information
	void *VP; // pointers to the no holders of other information
} mtype; // molecule types, with other info 
typedef struct { // types -- container for other type info
	int na; // # of atom types
	atype *a; // na of these
	int nr; // # of residue types
	rtype *r; // nr of these
	int nm; // # of molecule types
	mtype *m; // nm of these
	int nVP ; // number void pointers
	void *VP; // the void pointers
} types; // superstructure for typing information

/***************** structure parameter_set ********************/
/* parameters for atoms, molecules, residues, etc. */
typedef struct {
	int nAT; // number of atom types represented in this set
	atype *AT; // array of atom type structures
	int nRT; // number of residue types represented in this set
	rtype *RT; // array of residue type structures
	int nMT; // number of molecule types represented in this set
	mtype *MT; // array of molecule type structures
	int nBT;
	bond_type *BT;
	int nHBT; // for hydrogen bonds
	bond_type *HBT;
	int nNBT; // for non-bonded interactions
	bond_type *NBT;
	int nANT;
	angle_type *ANT;
	int nTRT;
	torsion_type *TRT;
} parameter_set;
/** @}*/

#endif
