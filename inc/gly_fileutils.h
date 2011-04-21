/** \file  gly_fileutils.h
\brief header file that loads utilities related to file I/O
*/
#if !defined(GLYLIB_FILE_UTILITIES)
#define GLYLIB_FILE_UTILITIES

#include <sys/types.h> ///< Required for directory manipulations
#include <dirent.h> ///< Required for directory manipulations
#include <unistd.h> ///< Required to get the current working directory

/** \addtogroup FILE_UTILS 
 * @{*/

/****************  STRUCTURES *********************/
/** A structure to hold sets of keywords and values that have been,
for example, read from a file.  */
typedef struct {
        int n; ///< The number of keys and values in this set
        char **K; ///< An array of n keys
        char **V; ///< An array of n values
} gly_keysvals;

typedef struct {
	char *N; /// the name of the file -- use strdup to allocate/copy
	FILE *F; /// the file pointer
} fileset;
/**< A structure for ease in passing sets of file info between functions */

/// A structure for a function similar to the perl notion of "slurping" a file
typedef struct {
	int n; /// number of lines
	char **L; /// each line
} fileslurp; 

/** List of files and directories generated for each state.  */
//typedef struct{
        //char *N; ///< Name of this directory
        //int nF; ///< Number of plain files in the directory
        //char **F; ///< Names of the nF files
        //int *iF; ///< nF utility integers (mark for deletion or saving, etc.)
        //int nD; ///< Number of sub-directories
        //glylib_directory *D; ///< nD structures for the sub-directories
        //int *iD; ///< nD utility integers
//} glylib_directory;



/****************  FUNCTIONS *********************/

/** Opens a file with error checking. 
This function will open a file and then check to see if the file actually opened.  
If not, it will write a descriptive error message and exit(1) the program. 
Doing this saves many, many core dumps. */
FILE *myfopen(const char* Name,const char* option);///< Opens a file with error checking.
//char glyopendir(char *pathToDir);
/** Like myfopen, but the freopen version. */
FILE *myfreopen(const char* Name,const char*,FILE* option); ///< Like myfopen, but the freopen version.

/** Complains (EEK!!!) if a file read fails.  Prints the following: "

Unexpected surprise during read of file FileDesc.

Exiting.

"*/
void  read_eek(const char *surprise, const char *FileDesc); ///< Complains (EEK!!!) if a file read fails.

 /** Complains (EEK!!!) if a file read fails. Prints the following: "

" */
void read_neek(const char *surprise,int line,int field); ///< Complains (EEK!!!) if a file read fails.

/** Complains (EEK!!!) if a file read fails. */
void read_fneek(const char *surprise, int line, int field, const char *FILENAME); ///< Complains (EEK!!!) if a file read fails.

/** Returns the number of times the string str is found in file F */
int sscan_file(FILE* F,const char* str); ///< Scans a file for instances of a string.

/** Returns the number of times the character c is found in file F */
int cscan_file(FILE* F,char c); 

/** Slurps in the contents of an entire file or directory */
fileslurp slurp_file(fileset F); 
//fileslurp slurp_directory(char *pathToDir);

/** Determines the current working directory. */
const char *gly_get_current_working_directory(void);

/** Reads in keyword/value pairs from a file, drops them in a
structure designed for the purpose.

Some notes: 
	- If ignore_hash is 0, everything from a hash mark (#) to the end of the line is ignored.
	- Any line not containing the separator SEP is ignored.
	- Keywords and values may not be separated by a newline.
	- If newline_sep is set to zero, everything* from the start of the line to the first instance 
		of SEP will be regarded as the KEYWORD and everything* from the first instance of SEP
		to the end of the line will be regarded as the VALUE. *Intial and final whitespace will 
		be removed from both keyword and value, but all else, including internal whitespace,
		will be retained.
	- If newline_sep is not set to zero, neither keywords nor values may contain whitespace.
	- SEP should be passed as a string complete with a null terminator.
*/

gly_keysvals get_keysvals_from_slurp(fileslurp F, const char *SEP, int newline_sep, int ignore_hash);

/** This is like get_keysvals_from_slurp (see) except it uses a 
	file and not a fileslurp. 

In fact, this function merely slurps in your file and then calls get_keysvals_from_slurp 
*/
gly_keysvals get_keysvals_from_file(fileset F, const char *SEP, int newline_sep, int ignore_hash);
/*@}*/

/*
   Deallocation routines
*/
void deallocateFileslurp(fileslurp *F);

#endif
