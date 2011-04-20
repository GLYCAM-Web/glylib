#include "../inc/mylib.h"
#include "../inc/molecules.h"
#include "../inc/AMBER/amber.h"

/* \file deallocate_amber_structures.c 
\addtogroup MEMORY_MANAGEMENT
\brief   Dealloction routines for structures relevant to AMBER.

Begun on 20100304 by BLFoley 

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
/*
void deallocateX(X *a){
 int i;
 //set null
 a[0].t = 0x0;
 // deallocate each
 if(a[0].b != NULL && a[0].b != 0x0){
 	for(i=0;i<a[0].nb;i++){deallocateX(&a[0].b[i]);}
 	free(a[0].b);}
 // free each
 if(a[0].OD != NULL && a[0].OD != 0x0){
 	for(i=0;i<a[0].nOD;i++){if(a[0].OD[i] != NULL && a[0].OD[i] != 0x0){free(a[0].OD[i]);}}
	free(a[0].OD);}
 // free top
 if(a[0].N != NULL && a[0].N != 0x0){free(a[0].N);}
 // check and warn if non-null
 if(a[0].VP!= NULL && a[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating X structure with non-NULL void pointer!\n");}
 return ;
}
*/

/********* structures from  AMBER/amber_prmtop_structs.h *****************/
void deallocateAmberPrmtopSection(amber_prmtop_section *aps){
 int i;
 // free each
 if(aps[0].D != NULL && aps[0].D != 0x0){
 	for(i=0;i<aps[0].nt;i++){if(aps[0].D[i] != NULL && aps[0].D[i] != 0x0){free(aps[0].D[i]);}}
	free(aps[0].D);}
 // free top
//printf("aps[0].N is %d\n",aps[0].N);
 if(aps[0].N != NULL && aps[0].N != 0x0){free(aps[0].N);}
 if(aps[0].FORMAT != NULL && aps[0].FORMAT != 0x0){free(aps[0].FORMAT);}
//printf("aps[0].TYPE is %d\n",aps[0].TYPE);
 if(aps[0].TYPE != NULL && aps[0].TYPE != 0x0){free(aps[0].TYPE);}
 if(aps[0].desc != NULL && aps[0].desc != 0x0){free(aps[0].desc);}
 return ;
}
void deallocateAmberPrmtop(amber_prmtop *ap){
 int i;
 // deallocate each
 if(ap[0].S != NULL && ap[0].S != 0x0){
 	for(i=0;i<ap[0].nS;i++){
//printf("Deallocating prmtop section i=%d; it's FLAG is >>>%s<<<\n",i,ap[0].S[i].N);
		deallocateAmberPrmtopSection(&ap[0].S[i]);}
 	free(ap[0].S);}
 // free each
 if(ap[0].SN != NULL && ap[0].SN != 0x0){
 	for(i=0;i<ap[0].nS;i++){if(ap[0].SN[i] != NULL && ap[0].SN[i] != 0x0){free(ap[0].SN[i]);}}
	free(ap[0].SN);}
 if(ap[0].FLAGS != NULL && ap[0].FLAGS != 0x0){
 	for(i=0;i<ap[0].nFF;i++){if(ap[0].FLAGS[i] != NULL && ap[0].FLAGS[i] != 0x0){free(ap[0].FLAGS[i]);}}
	free(ap[0].FLAGS);}
 // free top
 if(ap[0].FN != NULL && ap[0].FN != 0x0){free(ap[0].FN);}
 if(ap[0].VERSION != NULL && ap[0].VERSION != 0x0){free(ap[0].VERSION);}
 if(ap[0].SS != NULL && ap[0].SS != 0x0){free(ap[0].SS);}
 if(ap[0].ES != NULL && ap[0].ES != 0x0){free(ap[0].ES);}
 if(ap[0].FORMATS != NULL && ap[0].FORMATS != 0x0){free(ap[0].FORMATS);}
 return ;
}

/********* structures from  AMBER/read_prep.h *****************/
/*  Not writing deallocations for these right now. 
 *  No time... 20100304 BLFoley
 */
/*
typedef struct{
	char *FN; // file name from which this set of preps originated
	char *NAMDBF; // name of database to which this set belongs
	int IDBGEN , IREST , ITYPF; // Amber flags from top of file
} ambermprepinfo;
typedef struct{
	char *TITLE; // description of residue from prep file
	char *NAMF; // NAMF from residue
	char *NAMRES; // name of the residue from the prep file
	char *INTX; // coordinates are internal "INT" or xyz "XYZ"
	int  KFORM; // KFORM from amber file
	char *IFIXC; // geometry format
	char *IOMIT; // flag for dummy atom omission
	char *ISYMDU; // symbol used for dummy atoms
	char *IPOS; // flag for which dummy atoms to delete
	double CUT; // cutoff for assuming atoms are bonded
	int nIOPR; // number of IOPR cards for this residue
	llcharset *IOPR; // content of those cards ("DONE" is not included)
} amberrprepinfo;
typedef struct{
	int I;
	char *IGRAPH;
	char *ISYMBL;
	char *ITREE;
	int NA;
	int NB;
	int NC;
	double R;
	double THETA;
	double PHI;
	double CHG;
} amberaprepinfo;
*/

/********* structures from  AMBER/run_control.h *****************/
/*  Not writing deallocations for these right now. 
 *  No time... 20100304 BLFoley
 */
/*
typedef struct {
	int WANT_MDEN; ///< 0=YES ; 1=NO (make mden output file?)
	int MDEN_EVERY; ///< Step-frequency with which to write entries to mden file
	int WANT_AUTOCORR; ///< Should the program call get_autocorr_est_array()?
	char *AUTOCORR_DATA; ///< List of mden entries to get autocorrelations for
	int WANT_AUTOCORR_GRAPH; ///< 0=YES ; 1=NO (might require another program)
	int WANT_GRAPH; ///< 0=YES ; 1=NO (might require another program)
	char *GRAPH_DATA; ///< List of graphs to make -- all versus step number
} AMBER_MD_DATA_REQUEST;
typedef struct {
        int NUMBER; ///< A user-defined number; not the order found in the file
        char *COMMENT; ///< Free-form user comment about this change
        char *CHANGE_BY_WHICH;///< How to make the change, e.g.,RESIDUE_TYPE NOT-SOLVENT
        char *SOLVENT_ID;///< List of residues to consider as solvent
        char *RESIDUE_LIST; ///< List of residues to consider
        char *ATOM_LIST; ///< List of atoms to consider
        char *CHANGE_TYPE;///< how to change the numbers (ABSOLUTE, FRACTION, etc.)
        double CHANGE_QUANT;///< Meaning will depend on the value in CHANGE_TYPE
} EE_TI_CHANGE_INFO;
typedef struct {
        double KLAMBDA; ///< The value of klambda.  Typically integral, but might float...  
        char *CLAMBDA_HOW; ///< Method for assigning CLAMBDA's.  Not necessarily the same as the integration technique (see ANALYZE_HOW)
        //## CLAMBDA_HOW = GAUSSIAN 12
        //## CLAMBDA_HOW = SECTIONS 20 
        char *POINTSFILE; ///< If CLAMBDA_HOW involves reading a file, this is the file to read.
        int nCLAMBDA; ///< The number of CLAMBDA values.
        double *CLAMBDA;///< The CLAMBDA values.
        int nWIDTHS; ///< Should equal nCLAMBDA or zero -- here to facilitate checks to see if widths are defined.
        double *WIDTHS; ///< Widths to assign to each CLAMBDA value (e.g., for GAUSSIAN quadrature)
        char *ANALYZE_HOW; ///< The numerical integration method to use.
        //## ANALYZE_HOW = GAUSSIAN
        //## ANALYZE_HOW = TRAPEZOID
        //## ANALYZE_HOW = SIMPSONS_RULE
} TI_INTEGRATION_INFO;
typedef struct{ 
	char AO; ///< Append/Overwrite files ('A', 'O' or '\0')
	char *mdin; ///< input control data for the min/md run
	char *mdout; ///< output user readable state info and diagnostics 
	char *prmtop; ///< input molecular topology, force ﬁeld, periodic box type, atom and residue names
	char *inpcrd; ///< input initial coordinates and (optionally) velocities and periodic box size
	char *restrt; ///< output ﬁnal coordinates, velocity, and box dimensions if any - for restarting run
	char *refc; ///< input (optional) reference coords for position restraints; also used for targeted MD
	char *mdcrd; ///< output coordinate sets saved over trajectory
	char *inptraj; ///< input input coordinate sets in trajectory format, when imin=5
	char *mdvel; ///< output velocity sets saved over trajectory
	char *mden; ///< output extensive energy data over trajectory
	char *mdinfo; ///< output latest mdout-format energy info
	char *radii; ///< Obviously, some sort of radii, but not documented at this point in the manual... (BLF)
	char *inpdip; ///< input polarizable dipole ﬁle, when indmeth=3
	char *rstdip; ///< output polarizable dipole ﬁle, when indmeth=3
	char *cpin; ///< input protonation state deﬁnitions
	char *cprestrt; ///< protonation state deﬁnitions, ﬁnal protonation states for restart (same format as cpin)
	char *cpout; ///< output protonation state data saved over trajectory
	char *evbin; ///< input input for EVB potentials 
} sander_pmemd_command_line;
typedef struct {
	int nclambda; ///< The number of clambda values held in this struct
	double *clambda; ///< clambda values (see icfe and klambda) -- an array here b/c they are usually in groups
	double comp; ///< The compressibility of a system (default = 44.6);
	double cut; ///< Nonbonded cutoff (default = 8.0);
	double dielc; ///< Dielectric for electrostatics (not same as GB simulations) 
	double drms; ///< convergence criterion for minimization (see manual)
	double dt; ///< Time step (psec) for MD runs (recommended max is 0.0020);
	double dx0; ///< Initial min step length to use (program will self-adjust)
	int ibelly; ///< Belly type dynamics.  Default is 0 (don't use).  See manual.
	int icfe; ///< Flag for thermodynamic integration; no=0; yes=1; (see clambda and klambda)
	int igb; ///< Flag for using GB or PB IS models; See Manual.
	int imin; ///< no_min=0; min=1; read_traj_and_analyze=5
	int ioutfm; ///< crd & vel output formats; normal=0; NetCDF=1
	int irest; ///< don't_restart=0; restart=1
	int iwrap; ///< don't_wrap=0; wrap=1
	double klambda; ///< klambda value (see icfe and clambda)
	int maxcyc; ///< Maximum number of minimization cycles
	int ncyc; ///< at ncyc, switch from steepest descent to conjugate graduent.
	int nmropt; ///< none=0; some>0; NOESY=2
	int nsnb; ///< Frequencey of nonbonded list updates (default = 25);
	int nstlim; ///< Number of MD simulation steps
	int ntave; ///< Every ntave steps print averages over last ntave steps;  =0 to disable averaging
	int ntb; ///< =1 for const. volume ; =2 for constant pressure ; =0 for no periodicity 
	int ntc; ///< No Shake = 1; SHAKE=2; constrain all bonds = 3
	int ntf; ///< Set equal to ntc or see manual for more information
	int ntmin; ///< Min Method. See manual.  If =1, need to set ncyc for switch
	int ntp; ///< Sets pressure control
	int ntpr; ///< Every ntpr steps, print info to mdinfo and mdout
	int ntr; ///< Sets restraints information; =0 for none; >0 for some (see restraint* & refc)
	int ntrx; ///< refc format text=1; binary=0
	int ntt; ///< Sets temperature control Const.Energy=0; >0 is some other control (=1 usually, see manual)
	int ntwe; ///< Eery ntwe steps write compact energy to mden file; =0 to inhibit
	int ntwprt; ///< Only print atoms 1 through ntwprt in crd files; =0 to print all
	int ntwr; ///< Every ntwr steps write a restart file
	int ntwv; ///< Every ntwx velocities are written to mdvel; =0 to inhibit; =-1 to combine with crd file
	int ntwx; ///< Every ntwx coordinates are written to mdcrd; =0 to inhibit
	int ntx; ///< only_crd=1; crd+vel=5
	int ntxo; ///< format of restrt file (output); =1 only (=0 is deprecated)
	double pres0; ///< Desired constant pressure (bars);  default = 1.0; 
	double presi; ///< Initial pressure 
	double restraint_wt; ///< Sixe of k, the restraining spring constant.
	char *restraintmask; ///< List of atoms to be restrained if ntr=1
	double scee; ///< Divide 1-4 electrostatic interactions by scee ; GLYCAM needs = 1.0;
	double scnb; ///< Divide 1-4 vdW interactions by scnb; GLYCAm needs = 1.0;
	double taup; ///< Pressure relaxation time = 0.1;
	double tautp; ///< Time constant for heat-bath coupling; default = 1.0; 
	double temp0; ///< Desired constant temperature (kelvin) 
	double tempi; ///< Initial temperature
} sander_pmemd_cntrl;
typedef struct {
	int nbflag; ///< Tells how to determine when to update the nonbonded list.  See manual.
	double skinnb; ///< Used when nbflag=1.  See manual.
} sander_pmemd_ewald;
typedef struct {
	char *type; ///< the type of change, e.g., TEMP0, TAUTP, etc.
	int istep1,istep2; ///< beginning and ending steps for the change
	int value1,value2; ///< change from value1 at istep1 to value2 at istep2
} sander_pmemd_nmropt_wt_change;
typedef struct {
	char *Title; ///< The title of the run
	sander_pmemd_cntrl SPC; ///< Control list for the run
	int use_SPE; ///< include ewald info if 0; don't if 1
	sander_pmemd_ewald SPE; ///< 
	int nSPNWT; ///< number of wt change info sets to include (none if 0)
	sander_pmemd_nmropt_wt_change SPNWT; ///< entries in the "&wt" section for nmropt
} sander_pmemd_MD_control_info;
typedef struct{
	char *Title; ///< free-form descriptor
	char *exe; ///< executable name (with path, etc.)
	sander_pmemd_MD_control_info MDCI; ///< Control info for this run
	int nopt; ///< Number of command-line options
	char **optd; ///< The nopt dash options ('\0' if there is an argument in the next list without one of these)
	char **optt; ///< The nopt option strings ('\0' if there is an argument in the previous list without one of these)
	int nouttst; ///< Number of output tests to perform
	char **outdone; ///< String to check if job is complete
	char **outfail; ///< String to check if job has failed
	char **outdoing; ///< String to check if job might be still running and ok
	char **outpass; ///< String to check if job is finished and is probably "successful"
	char **outdesc; ///< Description of overall result from tests
	int **state; ///< Integer description of result: failed=-1; success=0; unknown or in progress=+1;
	int SUMMARY; ///< Integer summary: one or more failed=-1; all success=0; unknown or in progress=+1;
} sander_pmemd_MD_run_info;
*/
