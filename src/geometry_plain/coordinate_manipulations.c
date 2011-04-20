#include <mylib.h>
#include <molecules.h>
#include <geometries.h>

/** Euclidean distance from x to y
 */
double
get_distance_AB_points(coord_3D a, coord_3D b)
{ return sqrt(pow(b.i-a.i,2) + pow(b.j-a.j,2) + pow(b.k-a.k,2)); }


double get_angle_between_vectors(vectormag_3D a, vectormag_3D b)
{
	return acos(get_dotprod(normalize_vec(a),normalize_vec(b)));

}
double get_angle_ABC_points(coord_3D a, coord_3D b, coord_3D c)
{
	vectormag_3D ab = coord_to_vec(subtract_coord(a,b));
	vectormag_3D cb = coord_to_vec(subtract_coord(c,b));
	return get_angle_between_vectors(ab,cb);
}

/** This assumes that the bonding goes as: A-B-C-D or that the input vectors
	follow an equivalent direction convention.  When A and D are eclipsed, 
	the angle is zero.  */
double get_dihedral_ABCD_points(coord_3D a, coord_3D b, coord_3D c, coord_3D d)
{
	vectormag_3D b1 = coord_to_vec(subtract_coord(a, b));
	vectormag_3D b2 = coord_to_vec(subtract_coord(b, c));
	vectormag_3D b3 = coord_to_vec(subtract_coord(c, d));

	return atan2( get_dotprod(scalarmult_vec(b1, b2.d), get_crossprod(b2,b3)),
		get_dotprod(get_crossprod(b1,b2), get_crossprod(b2,b3)) ) ;
}


void
translate_coords_dp_list(coord_3D **coords, int num_coords, coord_3D shift)
{
    int i;

    for(i=0; i!=num_coords; i++){
        coords[i]->i += shift.i;
        coords[i]->j += shift.j;
        coords[i]->k += shift.k;
    }
}


void
orient_coords2_to_coords1_dp_list(coord_3D **coords, int num_coords,
          const coord_3D *bond_atom_a, coord_3D *bond_atom_b, double distance,
          const coord_3D *angle_atom_a, double theta,
          coord_3D *angle_atom_b, double rho,
          const coord_3D *dih_atom_a, coord_3D *dih_atom_b, double tau,
          const coord_3D *tor_atom_a, const coord_3D *ref_angle_a, double phi,
          coord_3D *tor_atom_b, coord_3D *ref_atom_b, double omega)
{
    coord_3D d = get_cartesian_point_from_internal_coords(*ref_angle_a, *tor_atom_a, *bond_atom_a,
                                    theta, PI - phi, distance);

    translate_coords_dp_list(coords, num_coords, subtract_coord(d,*bond_atom_b));

    double cur_dihedral = get_dihedral_ABCD_points(*dih_atom_a, *bond_atom_a,
                                       *bond_atom_b, *dih_atom_b);
    double cur_bond_angle = get_angle_ABC_points(*bond_atom_a, *bond_atom_b, *angle_atom_b);

    rotate_coords_about_axis_dp_list(coords, num_coords, *bond_atom_b,
                  get_crossprod(
                    get_vector_from_coords(*bond_atom_a, *bond_atom_b),
                    get_vector_from_coords(*angle_atom_b, *bond_atom_b)
                  ),
                  rho-cur_bond_angle);

    /*rotate_coords_about_axis_dp_list(coords, num_coords, *bond_atom_a,
                  subtract(*bond_atom_a, *bond_atom_b), cur_dihedral-tau);*/
    rotate_coords_about_axis_dp_list(coords, num_coords, *bond_atom_a,
                  get_vector_from_coords(*bond_atom_b, *bond_atom_a), cur_dihedral-tau);

    double cur_omega = get_dihedral_ABCD_points(*bond_atom_a, *bond_atom_b,
                                    *tor_atom_b, *ref_atom_b);

    /*rotate_coords_about_axis_dp_list(coords, num_coords, *bond_atom_b,
                  subtract(*bond_atom_b, *tor_atom_b),
                  cur_omega-omega);*/
    rotate_coords_about_axis_dp_list(coords, num_coords, *bond_atom_b,
                  get_vector_from_coords(*bond_atom_a, *tor_atom_b),
                  cur_omega-omega);
}

