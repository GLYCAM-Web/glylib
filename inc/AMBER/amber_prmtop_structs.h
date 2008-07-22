/** \file amber_prmtop_structs.h Structure for amber prmtop file info 

begun on 20080424 by BLFoley

If you add new flags, also add formats, preserving order -AND- define the
new flag in the amber_prmtop structure below and in the initialization
function (currently amber_prmtop_init.c).

NOTE: one space to separate fields, no space at end of set!
*/

#if !defined(GLYLIB_AMBER_PRMTOP_STRUCTS)
#define GLYLIB_AMBER_PRMTOP_STRUCTS

#define AMBER_PRMTOP_FLAGS "TITLE POINTERS ATOM_NAME CHARGE MASS \
ATOM_TYPE_INDEX NUMBER_EXCLUDED_ATOMS NONBONDED_PARM_INDEX RESIDUE_LABEL RESIDUE_ID \
RESIDUE_POINTER BOND_FORCE_CONSTANT BOND_EQUIL_VALUE ANGLE_FORCE_CONSTANT ANGLE_EQUIL_VALUE \
DIHEDRAL_FORCE_CONSTANT DIHEDRAL_PERIODICITY DIHEDRAL_PHASE SOLTY LENNARD_JONES_ACOEF \
LENNARD_JONES_BCOEF BONDS_INC_HYDROGEN BONDS_WITHOUT_HYDROGEN ANGLES_INC_HYDROGEN ANGLES_WITHOUT_HYDROGEN \
DIHEDRALS_INC_HYDROGEN DIHEDRALS_WITHOUT_HYDROGEN EXCLUDED_ATOMS_LIST HBOND_ACOEF HBOND_BCOEF \
HBCUT AMBER_ATOM_TYPE TREE_CHAIN_CLASSIFICATION JOIN_ARRAY IROTAT \
SOLVENT_POINTERS ATOMS_PER_MOLECULE BOX_DIMENSIONS CAP_INFO CAP_INFO2 \
RADIUS_SET RADII SCREEN LES_NTYP LES_TYPE \
LES_FAC LES_CNUM LES_ID PERT_BOND_ATOMS PERT_BOND_PARAMS \
PERT_ANGLE_ATOMS PERT_ANGLE_PARAMS PERT_DIHEDRAL_ATOMS PERT_DIHEDRAL_PARAMS PERT_RESIDUE_NAME \
PERT_ATOM_NAME PERT_ATOM_SYMBOL ALMPER IAPER PERT_ATOM_TYPE_INDEX \
PERT_CHARGE POLARIZABILITY PERT_POLARIZABILITY"
#define AMBER_PRMTOP_FORMATS "20a4 10I8 20a4 5E16.8 5E16.8 \
10I8 10I8 10I8 20a4 10a8 \
10I8 5E16.8 5E16.8 5E16.8 5E16.8 \
5E16.8 5E16.8 5E16.8 5E16.8 5E16.8 \
5E16.8 10I8 10I8 10I8 10I8 \
10I8 10I8 10I8 5E16.8 5E16.8 \
5E16.8 20a4 20a4 10I8 10I8 \
3I8 10I8 5E16.8 10I8 5E16.8 \
1a80 5E16.8 5E16.8 10I8 10I8 \
5E16 10I8 10I8 10I8 10I8 \
10I8 10I8 10I8 10I8 20a4 \
20a4 20a4 5E16.8 10I8 10I8 \
5E16.8 5E16.8 5E16.8"

// structure representing a section in a topology file
typedef struct {
	char *N; // the "FLAG" name of this section
	char *FORMAT; // the formatting information, as found in the file
	char *TYPE; // the type of data contained in this section, e.g., char, int, double (no float...)
	char *desc; // description of data in this section (free format)
	int is_standard; // =0 if a standard section, not =0 if not standard
	int nc; // number of characters per datum
	int npl; // data points per line
	int nt; // total number of pieces of data
	int ns; // number of datum points per set (e.g., =3 for x,y,z data)
	int nu; // number of units of information (=nt/ns) (one of these will get deleted... maybe)
	char **D; // pointer to nt bits of data (will read in the prmtop file as-is and convert later)
} amber_prmtop_section;

// structure amber_prmtop
//
// IMPORTANT NOTES 
// After the first few entries about filenames and headers, all the entries,
//	despite being named for AMBER prmtop variables, are NOT the actual
//	variables.  Rather, they are pointers to positions in the section
//	structures.
// Programmers wishing to alter/update this should also consider the contents
//	of the initialization function (currently amber_prmtop_init.c) to 
//	learn the data types used in the sections.
typedef struct {
	char *FN; // name of file from which the data were obtained
        char *VERSION; // the top line that usually discusses version info
	int nS; // number of sections found
	char **SN; // names of sections, in the order found
	amber_prmtop_section *S; // addresses for the actual sections
	int nSS; // number of standard sections found
	int *SS; // pointers to those standard sections
	int nES;// number of extra sections found
	int *ES; // list of ES addresses
	int nFF,*FORMATS; // for programming later
	char **FLAGS; // list of known flags
//FORMAT(20a4)  (ITITL(i), i=1,20)
        int ITITL; // the title section
// Read the "pointers" section in here, as well as into an 
//	amber_prmtop_section, but save room for extras if found
	int POINTERS; // pointer to the section containing the original char-strings
//FORMAT(12i6) 
	int NATOM  ; // total number of atoms 
  	int NTYPES ; // total number of distinct atom types
  	int NBONH  ; // number of bonds containing hydrogen
  	int MBONA  ; // number of bonds not containing hydrogen
  	int NTHETH ; // number of angles containing hydrogen
  	int MTHETA ; // number of angles not containing hydrogen
  	int NPHIH  ; // number of dihedrals containing hydrogen
  	int MPHIA  ; // number of dihedrals not containing hydrogen
  	int NHPARM ; // currently not used
  	int NPARM  ; // currently not used
  	int NEXT   ; // number of excluded atoms
  	int NRES   ; // number of residues
  	int NBONA  ; // MBONA + number of constraint bonds
  	int NTHETA ; // MTHETA + number of constraint angles
  	int NPHIA  ; // MPHIA + number of constraint dihedrals
  	int NUMBND ; // number of unique bond types
  	int NUMANG ; // number of unique angle types
  	int NPTRA  ; // number of unique dihedral types
  	int NATYP  ; // number of atom types in parameter file, see SOLTY below
  	int NPHB   ; // number of distinct 10-12 hydrogen bond pair types
  	int IFPERT ; // set to 1 if perturbation info is to be read in
  	int NBPER  ; // number of bonds to be perturbed
  	int NGPER  ; // number of angles to be perturbed
  	int NDPER  ; // number of dihedrals to be perturbed
  	int MBPER  ; // number of bonds with atoms completely in perturbed group
  	int MGPER  ; // number of angles with atoms completely in perturbed group
  	int MDPER  ; // number of dihedrals with atoms completely in perturbed groups
  	int IFBOX  ; // set to 1 if standard periodic box, 2 when truncated octahedral
  	int NMXRS  ; // number of atoms in the largest residue
  	int IFCAP  ; // set to 1 if the CAP option from edit was specified
// FORMAT(20a4)  (IGRAPH(i), i=1,NATOM)
  	int IGRAPH; // IGRAPH : the user atoms names 
// FORMAT(5E16.8)  (CHRG(i), i=1,NATOM)
  	int CHRG   ; // CHRG   : the atom charges.  
		// (Divide by 18.2223 to convert to charge in units of the electron charge) 
// FORMAT(5E16.8)  (AMASS(i), i=1,NATOM)
  	int AMASS  ; // AMASS  : the atom masses 
// FORMAT(12I6)  (IAC(i), i=1,NATOM)
  	int IAC    ; // IAC    : index for the atom types involved in Lennard Jones (6-12) 
           	// interactions.  See ICO below.  
// FORMAT(12I6)  (NUMEX(i), i=1,NATOM)
  	int NUMEX  ; // NUMEX  : total number of excluded atoms for atom "i".  See
           	// NATEX below.  
// FORMAT(12I6)  (ICO(i), i=1,NTYPES*NTYPES)
  	int ICO    ; // ICO    : provides the index to the nonbon parameter
           	// arrays CN1, CN2 and ASOL, BSOL.  All possible 6-12
           	// or 10-12 atoms type interactions are represented.
           	// NOTE: A particular atom type can have either a 10-12
           	// or a 6-12 interaction, but not both.  The index is
           	// calculated as follows:
             	// index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
           	// If index is positive, this is an index into the
           	// 6-12 parameter arrays (CN1 and CN2) otherwise it
           	// is an index into the 10-12 parameter arrays (ASOL and BSOL).  
// FORMAT(20A4)  (LABRES(i), i=1,NRES)
  	int LABRES ; // LABRES : the residue labels 
// ****** START HERE -- figure out the actual AMBER name for this...
// ------ LOOK vvvvvvvv
	int IRES;  // the residue numbers (of some sort)
// ------ LOOK ^^^^^^^^
// *****
// FORMAT(12I6)  (IPRES(i), i=1,NRES)
  	int IPRES  ; // IPRES  : atoms in each residue are listed for atom "i" in
           	// IPRES(i) to IPRES(i+1)-1 
// FORMAT(5E16.8)  (RK(i), i=1,NUMBND)
  	int RK     ; // RK     : force constant for the bonds of each type, kcal/mol 
// FORMAT(5E16.8)  (REQ(i), i=1,NUMBND)
  	int REQ    ; // REQ    : the equilibrium bond length for the bonds of each type, angstroms 
// FORMAT(5E16.8)  (TK(i), i=1,NUMANG)
  	int TK     ; // TK     : force constant for the angles of each type, kcal/mol A**2 
// FORMAT(5E16.8)  (TEQ(i), i=1,NUMANG)
  	int TEQ    ; // TEQ    : the equilibrium angle for the angles of each type, radians 
// FORMAT(5E16.8)  (PK(i), i=1,NPTRA)
  	int PK     ; // PK     : force constant for the dihedrals of each type, kcal/mol 
// FORMAT(5E16.8)  (PN(i), i=1,NPTRA)
  	int PN     ; // PN     : periodicity of the dihedral of a given type 
// FORMAT(5E16.8)  (PHASE(i), i=1,NPTRA)
  	int PHASE  ; // PHASE  : phase of the dihedral of a given type, radians 
// FORMAT(5E16.8)  (SOLTY(i), i=1,NATYP)
  	int SOLTY  ; // SOLTY  : currently unused (reserved for future use) 
// FORMAT(5E16.8)  (CN1(i), i=1,NTYPES*(NTYPES+1)/2)
  	int CN1    ; // CN1    : Lennard Jones r**12 terms for all possible atom type
           	// interactions, indexed by ICO and IAC; for atom i and j
           	// where i < j, the index into this array is as follows
           	// (assuming the value of ICO(index) is positive):
           	// CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).  
// FORMAT(5E16.8)  (CN2(i), i=1,NTYPES*(NTYPES+1)/2)
  	int CN2    ; // CN2    : Lennard Jones r**6 terms for all possible atom type
           	// interactions.  Indexed like CN1 above.  
/* NOTE: the atom numbers in the following arrays that describe bonds, 
	angles, and dihedrals are coordinate array indexes for runtime speed. 
	The true atom number equals the absolute value of the number divided by 
	three, plus one. In the case of the dihedrals, if the fourth atom is negative, 
	this implies that the dihedral is an improper. If the third atom is negative, 
	this implies that the end group interations are to be ignored. End group 
	interactions are ignored, for example, in dihedrals of various ring systems 
	(to prevent int counting of 1-4 interactions) and in multiterm dihedrals.  */
// FORMAT(12I6)  (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
  	int IBH    ; // IBH    : atom involved in bond "i", bond contains hydrogen
  	int JBH    ; // JBH    : atom involved in bond "i", bond contains hydrogen
  	int ICBH   ; // ICBH   : index into parameter arrays RK and REQ 
// FORMAT(12I6)  (IB(i),JB(i),ICB(i), i=1,NBONA)
  	int IB     ; // IB     : atom involved in bond "i", bond does not contain hydrogen
  	int JB     ; // JB     : atom involved in bond "i", bond does not contain hydrogen
  	int ICB    ; // ICB    : index into parameter arrays RK and REQ 
// FORMAT(12I6)  (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)
  	int ITH    ; // ITH    : atom involved in angle "i", angle contains hydrogen
  	int JTH    ; // JTH    : atom involved in angle "i", angle contains hydrogen
  	int KTH    ; // KTH    : atom involved in angle "i", angle contains hydrogen
  	int ICTH   ; // ICTH   : index into parameter arrays TK and TEQ for angle
           	// ITH(i)-JTH(i)-KTH(i) 
// FORMAT(12I6)  (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)
  	int IT     ; // IT     : atom involved in angle "i", angle does not contain hydrogen
  	int JT     ; // JT     : atom involved in angle "i", angle does not contain hydrogen
  	int KT     ; // KT     : atom involved in angle "i", angle does not contain hydrogen
  	int ICT    ; // ICT    : index into parameter arrays TK and TEQ for angle
           	// IT(i)-JT(i)-KT(i) 
// FORMAT(12I6)  (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i), i=1,NPHIH)
  	int IPH    ; // IPH    : atom involved in dihedral "i", dihedral contains hydrogen
  	int JPH    ; // JPH    : atom involved in dihedral "i", dihedral contains hydrogen
  	int KPH    ; // KPH    : atom involved in dihedral "i", dihedral contains hydrogen
  	int LPH    ; // LPH    : atom involved in dihedral "i", dihedral contains hydrogen
  	int ICPH   ; // ICPH   : index into parameter arrays PK, PN, and PHASE for
           	// dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i)

// FORMAT(12I6)  (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)
  	int IP     ; // IP     : atom involved in dihedral "i", dihedral does not contain hydrogen
  	int JP     ; // JP     : atom involved in dihedral "i", dihedral does not contain hydrogen
  	int KP     ; // KP     : atom involved in dihedral "i", dihedral does not contain hydrogen
  	int LP     ; // LP     : atom involved in dihedral "i", dihedral does not contain hydrogen
  	int ICP    ; // ICP    : index into parameter arrays PK, PN, and PHASE for
           	// dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i).  Note, if the
           	// periodicity is negative, this implies the following entry
           	// in the PK, PN, and PHASE arrays is another term in a
           	// multitermed dihedral.  
// FORMAT(12I6)  (NATEX(i), i=1,NEXT)
  	int NATEX  ; // NATEX  : the excluded atom list.  To get the excluded list for atom 
           	// "i" you need to traverse the NUMEX list, adding up all
           	// the previous NUMEX values, since NUMEX(i) holds the number
           	// of excluded atoms for atom "i", not the index into the 
           	// NATEX list.  Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
           	// excluded atoms are NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i)).  
// FORMAT(5E16.8)  (ASOL(i), i=1,NPHB)
  	int ASOL   ; // ASOL   : the value for the r**12 term for hydrogen bonds of all
           	// possible types.  Index into these arrays is equivalent
           	// to the CN1 and CN2 arrays, however the index is negative.
           	// For example, for atoms i and j, with i < j, the index is
           	// -ICO(NTYPES*(IAC(i)-1+IAC(j)).  
// FORMAT(5E16.8)  (BSOL(i), i=1,NPHB)
  	int BSOL   ; // BSOL   : the value for the r**10 term for hydrogen bonds of all
           	// possible types.  Indexed like ASOL.  
// FORMAT(5E16.8)  (HBCUT(i), i=1,NPHB)
  	int HBCUT  ; // HBCUT  : no longer in use 
// FORMAT(20A4)  (ISYMBL(i), i=1,NATOM)
  	int ISYMBL ; // ISYMBL : the AMBER atom types for each atom 
// FORMAT(20A4)  (ITREE(i), i=1,NATOM)
  	int ITREE  ; // ITREE  : the list of tree joining information, classified into five
           	// types.  M -- main chain, S -- side chain, B -- branch point, 
           	// 3 -- branch into three chains, E -- end of the chain 
// FORMAT(12I6)  (JOIN(i), i=1,NATOM)
  	int JOIN   ; // JOIN   : tree joining information, potentially used in ancient
           	// analysis programs.  Currently unused in sander or gibbs.  
// FORMAT(12I6)  (IROTAT(i), i = 1, NATOM)
  	int IROTAT ; // IROTAT : apparently the last atom that would move if atom i was
           	// rotated, however the meaning has been lost over time.
           	// Currently unused in sander or gibbs.
// The following are only present if IFBOX .gt. 0 
// FORMAT(12I6)  IPTRES, NSPM, NSPSOL
  	int IPTRES ; // IPTRES : final residue that is considered part of the solute,
           	// reset in sander and gibbs
  	int NSPM   ; // NSPM   : total number of molecules
  	int NSPSOL ; // NSPSOL : the first solvent "molecule" 
// FORMAT(12I6)  (NSP(i), i=1,NSPM)
  	int NSP    ; // NSP    : the total number of atoms in each molecule,
           	// necessary to correctly perform the pressure scaling.  
// FORMAT(5E16.8)  BETA, BOX(1), BOX(2), BOX(3)
  	int BETA   ; // BETA   : periodic box, angle between the XY and YZ planes in degrees.
  	int BOX    ; // BOX    : the periodic box lengths in the X, Y, and Z directions 
// The following are only present if IFCAP .gt. 0 
// FORMAT(12I6)  NATCAP
  	int NATCAP ; // NATCAP : last atom before the start of the cap of waters placed by edit 
// FORMAT(5E16.8)  CUTCAP, XCAP, YCAP, ZCAP
  	int CUTCAP ; // CUTCAP : the distance from the center of the cap to the outside
  	int XCAP   ; // XCAP   : X coordinate for the center of the cap
  	int YCAP   ; // YCAP   : Y coordinate for the center of the cap
  	int ZCAP   ; // ZCAP   : Z coordinate for the center of the cap 
// The following is only present if IFPERT .gt. 0
/* Note that the initial state, or equivalently the prep/link/edit state, 
	is represented by lambda=1 and the perturbed state, or final state 
	specified in parm, is the lambda=0 state. */ 
// FORMAT(12I6)  (IBPER(i), JBPER(i), i=1,NBPER)
  	int IBPER  ; // IBPER  : atoms involved in perturbed bonds
  	int JBPER  ; // JBPER  : atoms involved in perturbed bonds 
// FORMAT(12I6)  (ICBPER(i), i=1,2*NBPER)
  	int ICBPER ; // ICBPER : pointer into the bond parameter arrays RK and REQ for the
           	// perturbed bonds.  ICBPER(i) represents lambda=1 and 
           	// ICBPER(i+NBPER) represents lambda=0.  
// FORMAT(12I6)  (ITPER(i), JTPER(i), KTPER(i), i=1,NGPER)
  	int IPTER  ; // IPTER  : atoms involved in perturbed angles
  	int JTPER  ; // JTPER  : atoms involved in perturbed angles
  	int KTPER  ; // KTPER  : atoms involved in perturbed angles 
// FORMAT(12I6)  (ICTPER(i), i=1,2*NGPER)
  	int ICTPER ; // ICTPER : pointer into the angle parameter arrays TK and TEQ for 
           	// the perturbed angles.  ICTPER(i) represents lambda=0 and 
           	// ICTPER(i+NGPER) represents lambda=1.  
// FORMAT(12I6)  (IPPER(i), JPPER(i), KPPER(i), LPPER(i), i=1,NDPER)
  	int IPPER  ; // IPPER  : atoms involved in perturbed dihedrals
  	int JPPER  ; // JPPER  : atoms involved in perturbed dihedrals
  	int KPPER  ; // KPPER  : atoms involved in perturbed dihedrals
  	int LPPER  ; // LPPER  : atoms involved in pertrubed dihedrals 
// FORMAT(12I6)  (ICPPER(i), i=1,2*NDPER)
  	int ICPPER ; // ICPPER : pointer into the dihedral parameter arrays PK, PN and
           	// PHASE for the perturbed dihedrals.  ICPPER(i) represents 
           	// lambda=1 and ICPPER(i+NGPER) represents lambda=0.  
// FORMAT(20A4)  (LABRES(i), i=1,NRES)
  	int LRESPER ; // LABRES : residue names at lambda=0 
// FORMAT(20A4)  (IGRPER(i), i=1,NATOM)
  	int IGRPER ; // IGRPER : atomic names at lambda=0 
// FORMAT(20A4)  (ISMPER(i), i=1,NATOM)
  	int ISMPER ; // ISMPER : atomic symbols at lambda=0 
// FORMAT(5E16.8)  (ALMPER(i), i=1,NATOM)
  	int ALMPER ; // ALMPER : unused currently in gibbs 
// FORMAT(12I6)  (IAPER(i), i=1,NATOM)
  	int IAPER  ; // IAPER  : IAPER(i) = 1 if the atom is being perturbed 
// FORMAT(12I6)  (IACPER(i), i=1,NATOM)
  	int IACPER ; // IACPER : index for the atom types involved in Lennard Jones
           	// interactions at lambda=0.  Similar to IAC above.  See ICO above.  
// FORMAT(5E16.8)  (CGPER(i), i=1,NATOM)
  	int CGPER  ; // CGPER  : atomic charges at lambda=0 
// The following is only present if IPOL .eq. 1 
// FORMAT(5E18.8) (ATPOL(i), i=1,NATOM)
  	int ATPOL  ; // ATPOL  : atomic polarizabilities 
// The following is only present if IPOL .eq. 1 .and. IFPERT .eq. 1 
// FORMAT(5E18.8) (ATPOL1(i), i=1,NATOM)
  	int ATPOL1 ; // ATPOL1 : atomic polarizabilities at lambda = 1 (above is at lambda = 0) 
} amber_prmtop;

#endif
