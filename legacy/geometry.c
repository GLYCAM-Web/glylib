#include <stdlib.h>
#include "geometries.h"
//#include "geometry.h"

/* Magnitude
 */
/*double
norm(vector a)
{
    return sqrt(a.i*a.i + a.j*a.j + a.k*a.k);
}
*/

/*
vector
scale(vector v, double c)
{
    return (vector) {v.i*c,
                     v.j*c,
                     v.k*c};
}
*/

/* Cross product
 */
/*
vector
cross(vector a, vector b)
{
    return (vector) {a.j*b.k - a.k*b.j,
                     a.k*b.i - a.i*b.k,
                     a.i*b.j - a.j*b.i};
}
*/

/* Get the unit vector in given direction
 */
/*
vector
normalize(vector v)
{
    return scale(v, 1/norm(v));
}
*/

/* Dot product
 */
/*
double
dot(vector a, vector b)
{
    return a.i*b.i + a.j*b.j + a.k*b.k;
}
*/

/* Euclidean distance from x to y
 */
/*
double
dist(coord x, coord y)
{
    return sqrt(pow(y.i-x.i,2) + pow(y.j-x.j,2) + pow(y.k-x.k,2));
}
*/

/*
vector
subtract(vector a, vector b){
    return (vector) {a.i - b.i,
                     a.j - b.j,
                     a.k - b.k};
}
*/

/*
vector
add(vector a, vector b){
    return (vector) {a.i + b.i,
                     a.j + b.j,
                     a.k + b.k};
}
*/


double get_angle_two_vectors(vectormag_3D a, vectormag_3D b)
{
	return acos(get_dotprod(normalize_vec(a),normalize_vec(b)));

}
double get_angle_ABC(coord_3D a, coord_3D b, coord_3D c)
{
	vectormag_3D ab = coord_to_vec(subtract_coord(a,b));
	vectormag_3D cb = coord_to_vec(subtract_coord(c,b));
	return get_angle_between_vectors(ab,cb);
}

/*
double
get_angle(coord a, coord b, coord c)
{
    vector b1 = get_vector(b, a);
    vector b2 = get_vector(b, c);
    return acos( dot(b1,b2) / norm(b1) / norm(b2) );
}
*/


/** This assumes that the bonding goes as: A-B-C-D.  When A and D are eclipsed, 
	the angle is zero.  */
double get_dihedral_bonded_ABCD(coord_3D a, coord_3D b, coord_3D c, coord_3D d)
{
	vectormag_3D b1 = coord_to_vec(subtract_coord(a, b));
	vectormag_3D b2 = coord_to_vec(subtract_coord(b, c));
	vectormag_3D b3 = coord_to_vec(subtract_coord(c, d));

	return atan2( get_dotprod(scalarmult_vec(b1, b2.d), get_crossprod(b2,b3)),
		get_dotprod(get_crossprod(b1,b2), get_crossprod(b2,b3)) ) ;
}


/*
double
get_dihedral(coord a, coord b, coord c, coord d)
{
    vector b1 = get_vector(a, b);
    vector b2 = get_vector(b, c);
    vector b3 = get_vector(c, d);

    return atan2( dot(scale(b1, norm(b2)), cross(b2,b3)),
                  dot(cross(b1,b2), cross(b2,b3)) ) ;
}
*/


void
translate_coords_dp_list(coord **coords, int num_coords, coord shift)
{
    int i;

    for(i=0; i!=num_coords; i++){
        coords[i]->i += shift.i;
        coords[i]->j += shift.j;
        coords[i]->k += shift.k;
    }
}


void
orient_coords2_to_coords1_dp_list(coord **coords, int num_coords,
          const coord *bond_atom_a, coord *bond_atom_b, double distance,
          const coord *angle_atom_a, double theta,
          coord *angle_atom_b, double rho,
          const coord *dih_atom_a, coord *dih_atom_b, double tau,
          const coord *tor_atom_a, const coord *ref_angle_a, double phi,
          coord *tor_atom_b, coord *ref_atom_b, double omega)
{
    coord d = find_attachment_point(*ref_angle_a, *tor_atom_a, *bond_atom_a,
                                    theta, PI - phi, distance);

    translate_coords(coords, num_coords, subtract(d,*bond_atom_b));

    double cur_dihedral = get_dihedral(*dih_atom_a, *bond_atom_a,
                                       *bond_atom_b, *dih_atom_b);
    double cur_bond_angle = get_angle(*bond_atom_a, *bond_atom_b, *angle_atom_b);

    rotate_coords(coords, num_coords, *bond_atom_b,
                  cross(
                    get_vector(*bond_atom_a, *bond_atom_b),
                    get_vector(*angle_atom_b, *bond_atom_b)
                  ),
                  rho-cur_bond_angle);

    rotate_coords(coords, num_coords, *bond_atom_a,
                  subtract(*bond_atom_a, *bond_atom_b), cur_dihedral-tau);

    double cur_omega = get_dihedral(*bond_atom_a, *bond_atom_b,
                                    *tor_atom_b, *ref_atom_b);

    rotate_coords(coords, num_coords, *bond_atom_b,
                  subtract(*bond_atom_b, *tor_atom_b),
                  cur_omega-omega);
}

/*
coord **
atoms_to_coord_list(atom **atoms, int num_atoms)
{
    int i;

    coord **coords = (coord **)malloc(num_atoms * sizeof(coord *));
    for(i=0; i<num_atoms; i++)
        coords[i] = &atoms[i]->x;

    return coords;
}
*/
