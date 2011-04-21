/* \file geometries.h 
\brief Purpose: structures related to geometries of molecules,
	residues, atoms, etc.  Loaded when molecules.h is included.

begun by BLFoley on 20080606
   This file is loaded with the molecules.h header file by default.
        It should only be loaded explicitly in a program if, for some reason,
        you are using it separately. */
#if !defined(GLYLIB_GEOMETRIES)
#define GLYLIB_GEOMETRIES

/** \addtogroup  GEOMETRY
 * @{
 */
#if !defined(PI)
	#define PI 3.1415926535897932384626433832795028
#endif

#define get_vector_from_coords(a,b) coord_to_vec(subtract_coord(b,a))

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


/**************************************************************//**
			Functions
******************************************************************/
//void rotate_vector_to_Y_list(coord_3D*,int,vectormag_3D); // int is n
//void rotate_vector_to_X_list(coord_3D*,int,vectormag_3D); // int is n
//void rotate_vector_to_V_list(coord_3D*,int,vectormag_3D,vectormag_3D); // int is n 
void rotate_vector_to_Z_list(coord_3D*,int,vectormag_3D); // int is n
coord_3D get_geometric_center(coord_3D *c,int nc); // int is n
coord_3D get_geometric_center_dp(coord_3D **c,int nc); // int is n
plane get_plane(coord_3D,coord_3D,coord_3D);
plane get_plane_for_ring(int n,coord_3D *r); // gets the C-P average plane, atoms must be in order
vectormag_3D normalize_vec(vectormag_3D);
vectormag_3D scalarmult_vec(vectormag_3D,double);
vectormag_3D add_vec(vectormag_3D,vectormag_3D); // add two vectors
vectormag_3D subtract_vec(vectormag_3D,vectormag_3D); // subtract second vector from first
coord_3D scalarmult_coord(coord_3D,double); // add two vectors
coord_3D add_coord(coord_3D,coord_3D); // add two vectors
coord_3D subtract_coord(coord_3D,coord_3D); // subtract second vector from first
vectormag_3D get_crossprod(vectormag_3D, vectormag_3D); // returns the cross product
double get_dotprod(vectormag_3D, vectormag_3D); // returns dot product of two vectors
/* the following essentially returns the cosine of the angle between two vectors */
//double get_dotprodN(vectormag_3D, vectormag_3D); // returns dot prod, but normalizes vecs first
double get_magnitude(vectormag_3D); // calculates vector magnitude for "d" in structure
vectormag_3D zero_vec(); // zeros a vector
coord_3D zero_coord(); // zeros a coordinate set
coord_3D vec_to_coord(vectormag_3D); // turns a vector into a coordinate set
vectormag_3D coord_to_vec(coord_3D); // turns a coordinate set into a vector
void initialize_coord_3D(coord_3D *c);
void initialize_vectormag_3D(vectormag_3D *v);
void initialize_plane(plane *p);

/** Get angle between two vectors
 */
double get_angle_between_vectors(vectormag_3D a, vectormag_3D b);

/** Get angle abc
 */
double get_angle_ABC_points(coord_3D a, coord_3D b, coord_3D c);

/** Get the dihedral angle between planes abc and bcd
 */
double get_dihedral_ABCD_points(coord_3D a, coord_3D b, coord_3D c, coord_3D d);

/** Euclidean distance from x to y
 */
double get_distance_AB_points(coord_3D a, coord_3D b);

/** Translate a list of coordinates
 */
void translate_coords_dp_list(coord_3D **coords, int num_coords, coord_3D shift);

/** Create a rotation matrix for a rotation of theta degrees
 * about an axis through the given point in the given direction
 */
double *create_rotation_matrix(coord_3D point, vectormag_3D direction, 
                               double theta);
/** Deallocate the rotation matrix
 */
void destroy_rotation_matrix(double *matrix);

/** Apply the rotation matrix to the given coordinate
 */
void apply_rotation_matrix_to_coord_p(coord_3D *c, double *matrix);

/** Rotate a list of coordinates theta degrees about an axis through 
 * the given coordinate in the given direction
 */
void rotate_coords_about_axis_dp_list(coord_3D **coords, int num_coords, coord_3D point,
                   vectormag_3D direction, double theta);

/** Calculate coordinate d, where 
 * (1) d is distance units from c
 * (2) angle bcd is theta
 * (3) dihedral between planes abc and bcd is phi
 *
 * Link to derivation: 
 *
 * Note: should probably change the name of this
 */
coord_3D get_cartesian_point_from_internal_coords(coord_3D a, coord_3D b, coord_3D c,
                            double theta, double phi, double distance);

/** Note: change name/parameter names 
 */
void orient_coords2_to_coords1_dp_list(coord_3D **coords, int num_coords,
          const coord_3D *bond_atom_a, coord_3D *bond_atom_b, double distance,
          const coord_3D *angle_atom_a, double theta,
          coord_3D *angle_atom_b, double rho,
          const coord_3D *dih_atom_a, coord_3D *dih_atom_b, double tau,
          const coord_3D *tor_atom_a, const coord_3D *ref_angle_a, double phi,
          coord_3D *tor_atom_b, coord_3D *ref_atom_b, double omega);
/** @}*/

#endif
