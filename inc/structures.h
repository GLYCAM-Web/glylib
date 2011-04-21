/** \file structures.h 
\brief Very poorly named file containing structures
pertinent to the loading of pdb files.  Recommend name change.

Begun a long time ago by BLFoley and if it was modified at some point since 
then, the person who did so was likely Mike Tessier. */

#if !defined(GLYLIB_PDB_STRUCTURES)
#define GLYLIB_PDB_STRUCTURES

/** \addtogroup PDB
 * @{
 */
typedef struct {
	int a; /* integer describing the class for this record name */
	int b; /* integer describing record type in class */
} linetype;
typedef struct {
	char c[81]; /* pointers to characters in field */
} linefield;
typedef struct {
	char f[250]; // "for" string
	char r[250]; // "report" string
} queery_set;
typedef struct { 
	int a; // integer describing the class for this record name
	int b; // integer describing record type in class
	char m; /* Character to override modify flags */
	int n; /* the line number in the file */
	linefield f[25]; /* pointer to field */
	int ignore; /* if = 0 line already printed.  Ignore */
} linedef;
typedef struct {
        char c[50]; // for sets of short strings
} charlist;

typedef struct {
	char typ[100]; /* type of line -- padded for extra info.  Also
			used for line rule designations. */
	int f; /* number of fields */
	/* Since the length of any one line is 80 characters, there can be,
	*         at most (one hopes) 80 fields.  */
	int c[80]; /* number of characters in field */
	/*fieldinfo fi[80];  change information for each field */
	int i; /* toggle -- 0 for no change information, 1 if exists.
			Others designations can be added... */
	char t[80]; /* type of entry for each field:
	integer = i
	character = c
	string = s
	float = f
	double = d
	*/
} pdb_line_info;
typedef struct {
/* This structure defines a pdb line.  It defines it by class as given
  in the pdb format instructions.  Doing it this way allows rules to be
  set per class -- for example, one class should be present, at most, once
  in any pdb file.  */
	pdb_line_info b[25]; /* pointer to record type -- currently, there
			      are 22 record types in the most populated
			      class of records. */ 
	char r[50]; /* String containing a rule for this entire set of
		       lines. */
} pdb_line;

/* Begin declarations for amber line structure definitions */
pdb_line pdb_a[6]; /* There are 6 record classes, currently (20051104)*/
/* Begin declarations for modify-file rules.  The rules, for now, are limited
 * to rules involving HETATM entries.  Eventually, there will be more. */
pdb_line mod_a[6];
/* Pointer to a structure to keep track of user-specified modifications 
 * -- will be allocated in init_struct */
pdb_line umod_a[6];
/* Declare a pointer to a linedef.  This will be allocated in init_struct */
linedef *ln;
queery_set QS; //only one allowed currently
/** @}*/

#endif
