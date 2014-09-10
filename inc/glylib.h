#ifdef __cplusplus
extern "C" {
#endif
/** \file  glylib.h
\brief  A grand header file designed to simplify use of glylib by requiring 
that only one header file be loaded to get most of the glylib functionality.

	Begun on 20100216 by BLFoley

*/

#if !defined(GLYLIB_MAIN_HEADERS)
#define GLYLIB_MAIN_HEADERS

/** \addtogroup MAIN
 * @{
 */
/**************************************************************//**
	Standard GLYLIB structures and functions
******************************************************************/
#include <mylib.h> ///< Oft-used standard C headers and a few general ones from GLYLIB
#include <molecules.h> ///< Ensembles, molecules, etc.
#include <modes.h> ///< For vibrational modes (and similar)
#include <stats.h> ///< For calculating statistical information
#include <declarations.h> ///< Some declarations that don't fit easily elsewhere

/**************************************************************//**
	Structures and functions related to PDB files
******************************************************************/
#include <PDB.h>  /**< mostly for writing pdb's */
#include <load_pdb.h> ///< Specific to pdb files
    
/**************************************************************//**
	Structures and functions related to AMBER
******************************************************************/
#include <AMBER/amber.h> ///< Specific to AMBER
/** @}*/

#endif
#ifdef __cplusplus
}
#endif
