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
#if !defined(GLYLIB_MYLIB)
#define GLYLIB_MYLIB
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
#include <gly_fileutils.h>
#include <gly_codeutils.h>

/**
Defines that are are only needed for some specialized circumstances
*/
#define _GNU_SOURCE

/** Determines the current working directory. */
const char *gly_get_current_working_directory(void);

#endif
