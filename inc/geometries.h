/* File geometries.h begun by BLFoley on 20080606
   Purpose: structures related to geometries of molecules,
	residues, atoms, etc.  
   NOTE!!!  This file is loaded with the molecules.h header file by default.
        It should only be loaded explicitly in a program if, for some reason,
        you are using it separately. */
#if !defined(GLYLIB_GEOMETRIES)
#define GLYLIB_GEOMETRIES
typedef struct {
	double i,j,k; 
} coord_3D; // cartesian, polar, direction cosines, etc.

typedef struct {
	double i,j,k,d; // i,j,k vector with magnitude d
} vectormag_3D; // same as "plane" but with easier naming

typedef struct {
	double A,B,C,D;
} plane; // standard plane, Ax+By+Cz+D=0
#endif
