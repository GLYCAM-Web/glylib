#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>

typedef coord_3D vector;

typedef coord_3D coord;

/* Get angle abc
 */
double get_angle(coord a, coord b, coord c);

/* Get the dihedral angle between planes abc and bcd
 */
double get_dihedral(coord a, coord b, coord c, coord d);

/* Magnitude
 */
double norm(vector a);

/* Multiply vector by scalar
 */
vector scale(vector v, double c);

/* Cross product
 */
vector cross(vector a, vector b);

/* Get the unit vector in given direction
 */
vector normalize(vector v);

/* Dot product
 */
double dot(vector a, vector b);

/* Euclidean distance from x to y
 */
double dist(coord x, coord y);

vector subtract(vector a, vector b);

vector
add(vector a, vector b);

#define get_vector(a,b) subtract(b,a)

/* Translate a list of coordinates
 */
void translate_coords(coord **coords, int num_coords, coord shift);

/* Create a rotation matrix for a rotation of theta degrees
 * about an axis through the given point in the given direction
 */
double *create_rotation_matrix(coord point, vector direction, 
                               double theta);

/* Deallocate the rotation matrix
 */
void destroy_rotation_matrix(double *matrix);

/* Apply the rotation matrix to the given coordinate
 */
void apply_rotation_matrix(coord *c, double *matrix);

/* Rotate a list of coordinates theta degrees about an axis through 
 * the given coordinate in the given direction
 */
void rotate_coords(coord **coords, int num_coords, coord point,
                   vector direction, double theta);

/* Calculate coordinate d, where 
 * (1) d is distance units from c
 * (2) angle bcd is theta
 * (3) dihedral between planes abc and bcd is phi
 *
 * Link to derivation: 
 *
 * Note: should probably change the name of this
 */
coord find_attachment_point(coord a, coord b, coord c,
                            double theta, double phi, double distance);

/* Note: change name/parameter names 
 */
void fix_coords(coord **coords, int num_coords,
          const coord *bond_atom_a, coord *bond_atom_b, double distance,
          const coord *angle_atom_a, double theta,
          coord *angle_atom_b, double rho,
          const coord *dih_atom_a, coord *dih_atom_b, double tau,
          const coord *tor_atom_a, const coord *ref_angle_a, double phi,
          coord *tor_atom_b, coord *ref_atom_b, double omega);

/* Get a list of pointers to coordinates within the given list of atoms.
 * The list needs to be freed.
 */
//coord **atoms_to_coord_list(atom **atoms, int num_atoms);

#endif
