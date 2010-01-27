/* File geometries.h begun by BLFoley on 20080606
   Purpose: structures related to geometries of molecules,
	residues, atoms, etc.  
   NOTE!!!  This file is loaded with the molecules.h header file by default.
        It should only be loaded explicitly in a program if, for some reason,
        you are using it separately. */
#if !defined(GLYLIB_GEOMETRIES)
#define GLYLIB_GEOMETRIES

/**************************************************************//**
			Coordinates 
******************************************************************/
// Fixed Dimension
typedef struct {
	double i,j,k; ///< Values assigned to the three coordinates.
} coord_3D; ///< Any 3D coordinate: cartesian, polar, direction cosines, etc.  
typedef struct {
	double i,j,k,d; ///< i,j,k vector with magnitude d
} vectormag_3D; ///< Vector manipulations often require magnitudes, so this struct contains that information.

// Multi-Dimensional
typedef struct {
	int nD; ///< Number of dimensions.
	double *D;///< Values in each of the nD dimensions
} coord_nD; ///< An n-dimensional coordinate 
typedef struct {
	int nDi; ///< Number of dimensions.
	int *Di;///< Indicies for values in each of the nD dimensions
} nD_index; ///< Set of indices to n-dimensional coordinate
typedef struct {
	int nDp; ///< Number of dimensions.
	double **Dp;///< Pointers to values in each of the nD dimensions
} nD_ptrs; ///< Set of pointers to n-dimensional coordinate
typedef struct {
	int nD; ///< Number of dimensions.
	double *D;///< Values in each of the nD dimensions
	int *Di;///< Indicies for values in each of the nD dimensions
} coord_nDi; ///< An n-dimensional coordinate 



/**************************************************************//**
			Geometrical Objects 
******************************************************************/
typedef struct {
	double A,B,C,D;
} plane; ///< Standard cartesian plane with equation Ax+By+Cz+D=0

/**************************************************************//**
			Special Objects 
******************************************************************/
typedef struct { 
	char *STYPE,*GTYPE; ///< Types relevant to simulations (e.g., periodic) and to basic geometry (e.g., cubic)
	int nC; ///< Number of coordinates defined
	coord_nD *C; ///< nC coordinates (the dimensions are defined inside the structures, and might all be different)
	int nCD; ///< Number of descriptions defined
	char **CD; ///< nCD descriptions of the coordinate relevances (e.g, "lower corner" or "A in z=f(A)"). 
} boxinfo; ///< structure for holding (periodic or not) box information
typedef struct {
	double d; ///< a distance
	double a; ///< a relevant angle
	double t; ///< a reference torsion angle (sets chirality in the connection tree)
} bonded_position_set; ///< structure for info about a particular atom in terms of internal coordinates
#endif
