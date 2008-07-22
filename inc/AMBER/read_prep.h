/** \file read_prep.h Header file for function read_prep that reads the 
contents of a prep file into an existing molecule structure 

begun on 20071120 by BLFoley

These are loaded in the source
 - #include <mylib.h>  
 - #include <molecules.h>
 - #include <general.h>

See also AMBER documentation for prep file entry information
*/

#if !defined(GLYLIB_AMBER_PREP_HEADERS)
#define GLYLIB_AMBER_PREP_HEADERS

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


void read_prep(molecule *M, fileset F, types *TYP);
residue read_prepres(fileset F, types *TYP);
atom read_prepatom(const char *line, types *TYP);

#endif
