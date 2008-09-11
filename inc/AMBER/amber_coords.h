/** \file  AMBER/amber_coords.h
	\brief header file that loads utilities related to file I/O
*/
#if !defined(GLYLIB_AMBER_COORDS)
#define GLYLIB_AMBER_COORDS
/** Counts snapshots in a crd file

Makes certain that the number of snapshots read from a crd file makes
sense in terms of the given topology file.  This declaration should be moved to
some other location eventually (not so general -- really just amber). */
int   crdsnaps(const char *, const char *, int);
#endif
