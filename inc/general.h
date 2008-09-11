/** \file general.h 
\brief	Deprecated header file containing structures used by various functions
	within the library.  

This file now only calls the other headers to which its contents were moved.
*/

// Header file written by Lachele Foley, summer 2007

#if !defined(GLYLIB_GENERAL)
#define GLYLIB_GENERAL

#include <gly_codeutils.h>
#include <gly_fileutils.h>

/* An array of characters for times when a 2D char array might be confusing */
//typedef struct { char *T; } llcharset; //<<<< Moved to gly_codeutils.h

/* >>Convenient, but dreadfully wasteful, arrays of strings ("Text")
	>>Need get rid of these ASAP...  Can they go now?
Going now...  20080908 BLFoley  
DELETE THESE after, say, 20081201*/
//typedef struct { char T[501]; } lcharset;
//typedef struct { char T[251]; } mcharset;
//typedef struct { char T[51]; } scharset;
//typedef struct { char T[11]; } sscharset;

#endif
