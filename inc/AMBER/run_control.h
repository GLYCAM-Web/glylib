/** \file AMBER/run_control.h Contains information pertinent to amber simulations.

	Begun on 20080715 by BLFoley.  
 */

#if !defined(GLYLIB_AMBER_RUN_CONTROL)
#define GLYLIB_AMBER_RUN_CONTROL

/*! The following structure contains 
        lists of information requested by a user to be output from an amber
        simulation.  At present, however, it will be limited not only to 
        certain types of TI runs, but also to the specific way Lachele happens
        to want to run them.  If added to amber.h, it should be expanded.
 */
typedef struct {
	int WANT_MDEN; ///< 0=YES ; 1=NO (make mden output file?)
	int MDEN_EVERY; ///< Step-frequency with which to write entries to mden file
	int WANT_AUTOCORR; ///< Should the program call get_autocorr_est_array()?
	char *AUTOCORR_DATA; ///< List of mden entries to get autocorrelations for
	int WANT_AUTOCORR_GRAPH; ///< 0=YES ; 1=NO (might require another program)
	int WANT_GRAPH; ///< 0=YES ; 1=NO (might require another program)
	char *GRAPH_DATA; ///< List of graphs to make -- all versus step number
} AMBER_MD_DATA_REQUEST;


/*! This moves to amber.h eventually.  It will be used generally for TI calcs.
        Contains information on which charges to change and how to change them.
 */
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

/*! This will be used generally for TI calcs.
        Contains information about how the TI integration is to be carried out.
        Information includes beginning and end of run information.  For example,
        it contains the lambda values to use as well as the manner in which they
        should be integrated once the data are collected.
 */
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

/*! Contains the command-line requirements for a sander/pmemd run. 

sander [-help] [-O] [-A] -i mdin -o mdout -p prmtop -c inpcrd -r restrt
-ref refc -x mdcrd -y inptraj -v mdvel -e mden -inf mdinfo -radii radii
-cpin cpin -cpout cpout -cprestrt cprestrt -evbin evbin

-O Overwrite output files if they exist.
-A Append output files if they exist, (used mainly for replica exchange).
*/
typedef struct{ 
	char AO; ///< Append/Overwrite files ('A', 'O' or '\0')
	char *mdin; ///< input control data for the min/md run
	char *mdout; ///< output user readable state info and diagnostics 
		///< -o stdout will send output to stdout (to the terminal) instead of to a ﬁle.
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

/*! Contains "cntrl" information pertinent to a sander/pmemd run.  Currently,
	in the interest of time, I'm only writing it to contain the pieces I 
	typically use in a run.  The rest can be added later as these
	structs are easily extensible.
 */
typedef struct {
//
	int nclambda; ///< The number of clambda values held in this struct
	double *clambda; ///< clambda values (see icfe and klambda) -- an array here b/c they are usually in groups
	double comp; ///< The compressibility of a system (default = 44.6);
	double cut; ///< Nonbonded cutoff (default = 8.0);
//	
	double dielc; ///< Dielectric for electrostatics (not same as GB simulations) 
	double drms; ///< convergence criterion for minimization (see manual)
	double dt; ///< Time step (psec) for MD runs (recommended max is 0.0020);
	double dx0; ///< Initial min step length to use (program will self-adjust)
//
	int ibelly; ///< Belly type dynamics.  Default is 0 (don't use).  See manual.
	int icfe; ///< Flag for thermodynamic integration; no=0; yes=1; (see clambda and klambda)
	int igb; ///< Flag for using GB or PB IS models; See Manual.
	int imin; ///< no_min=0; min=1; read_traj_and_analyze=5
	int ioutfm; ///< crd & vel output formats; normal=0; NetCDF=1
	int irest; ///< don't_restart=0; restart=1
	int iwrap; ///< don't_wrap=0; wrap=1
//
	double klambda; ///< klambda value (see icfe and clambda)
//
	int maxcyc; ///< Maximum number of minimization cycles
//	
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
//
	double pres0; ///< Desired constant pressure (bars);  default = 1.0; 
	double presi; ///< Initial pressure 
//
	double restraint_wt; ///< Sixe of k, the restraining spring constant.
	char *restraintmask; ///< List of atoms to be restrained if ntr=1
//	
	double scee; ///< Divide 1-4 electrostatic interactions by scee ; GLYCAM needs = 1.0;
	double scnb; ///< Divide 1-4 vdW interactions by scnb; GLYCAm needs = 1.0;
//	
	double taup; ///< Pressure relaxation time = 0.1;
	double tautp; ///< Time constant for heat-bath coupling; default = 1.0; 
	double temp0; ///< Desired constant temperature (kelvin) 
	double tempi; ///< Initial temperature
//	
} sander_pmemd_cntrl;

/*! Contains control-ewald informatin pertinent to a sander/pmemd run.
	Currently, in the interest of time, I'm only writing it to contain 
	the pieces I typically use in a run.  The rest can be added later 
	as these structs are easily extensible. */
typedef struct {
	int nbflag; ///< Tells how to determine when to update the nonbonded list.  See manual.
	double skinnb; ///< Used when nbflag=1.  See manual.
} sander_pmemd_ewald;

/*! Contains control-change informatin pertinent to a sander/pmemd run.
	Currently, in the interest of time, I'm only writing it to contain 
	the pieces I typically use in a run.  The rest can be added later 
	as these structs are easily extensible. */
typedef struct {
	char *type; ///< the type of change, e.g., TEMP0, TAUTP, etc.
	int istep1,istep2; ///< beginning and ending steps for the change
	int value1,value2; ///< change from value1 at istep1 to value2 at istep2
} sander_pmemd_nmropt_wt_change;

/*! Contains information pertinent to a sander/pmemd run.  Currently, in the 
	interest of time, I'm only writing it to contain the pieces I 
	typically use in a run.  The rest can be added later as these
	structs are easily extensible.
 */
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


#endif
