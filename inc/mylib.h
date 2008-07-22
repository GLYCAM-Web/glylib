/** \file  mylib.h
	\brief generic header file that loads utilities 
		often needed by programmers

The header files included here are basic C library headers
(such as string.h). 
The other entries are declarations of functions that are expected
to be commonly used by programmers using the library.  For example,
functions for opening files and writing error messages are declared
in this file.  
*/
#if !defined(GLYLIB_HEAD)
#define GLYLIB_HEAD
/** 
Some libraries that are often needed 
*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <general.h>

/**
Defines that are often needed
*/
#define _GNU_SOURCE

/** 
Some general-use functions 
*/
/** Prints usage.*/
void  myusage(const char *USAGE); ///< Prints usage message and exits(1).
/** Prints error message and exits(1). 
Prints the following: 

"This program encountered the following error and will exit:

	whinetext " 
*/
void  mywhine(const char *whinetext); ///< Prints error message and exits(1).
/** Prints usage and an error message and exits(1). See mywhine for more information.
 Prints a usage statement after the whinetext*/
void  mywhineusage(const char *whinetext,const char *USAGE);///< Prints error message and usage and exits(1).
/** \defgroup FILE_UTILS "File Utilities" */
/** \addtogroup FILE_UTILS */
/*@{*/
/** Opens a file with error checking. 
This function will open a file and then check to see if the file actually opened.  
If not, it will write a descriptive error message and exit(1) the program. 
Doing this saves many, many core dumps. */
FILE *myfopen(const char* Name,const char* option);///< Opens a file with error checking.
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
/*! Makes certain that the number of snapshots read from a crd file makes
sense in terms of the given topology file.  This declaration should be moved to
some other location eventually (not so general -- really just amber). */
int   crdsnaps(const char *, const char *, int);///< counts snapshots in a crd file
/** 
* Usage:

*      const char * mylocaltimestring,mytime(format)

*              -for example:

*      strcpy(mylocaltimestring,mytime("format"));

* where:

*      format is one of:

*              [keyword]       [example]

*              YYYYMMDD        20070321

*              YYYYMMDDhhmmss  20070321152110

*              pretty          Wed, 21 March, 2007

*              longpretty      15:12:10, Wed, 21 March, 2007
*/
const char * mytime(const char*); ///< Gets the time in a nice format
/*! Returns the number of times the string str is found in file F */
int sscan_file(FILE* F,const char* str); ///< Scans a file for instances of a string.
/*! Returns the number of times the character c is found in file F */
int cscan_file(FILE* F,char c); ///< Scans a file for instances of a character.
fileslurp slurp_file(fileset F); ///< Slurps in the contents of an entire file
/*@{*/
#endif
