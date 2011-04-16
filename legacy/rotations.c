/* 
 * Functions for rotation about an axis
 *
 * author: Robert Davis
 * 	integrated with glylib by BLFoley 20110405
 */

#include <mylib.h>
#include <molecules.h>
#include <geometries.h>

/** Create a rotation matrix for a rotation of theta degrees
 * about an axis through the given point in the given direction
 */
double *
create_rotation_matrix(coord_3D point, vectormag_3D direction, double theta)
{
    direction=normalize_vec(direction);
    if(direction.d==0){fprintf(stderr,"Attempting to rotate about a vector of zero length in create_rotation_matrix.  Ignoring.\n");}
    double u = direction.i;
    double v = direction.j;
    double w = direction.k;

    double a = point.i;
    double b = point.j;
    double c = point.k;

    double *matrix = (double *)malloc(12 * sizeof(double));
    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;
    double cosT = cos(theta);
    double sinT = sin(theta);
 
    matrix[3] = a*(v2 + w2) - u*(b*v + c*w)
                 + (u*(b*v + c*w) - a*(v2 + w2))*cosT + (b*w - c*v)*sinT;

    matrix[7] = b*(u2 + w2) - v*(a*u + c*w)
                + (v*(a*u + c*w) - b*(u2 + w2))*cosT + (c*u - a*w)*sinT;

    matrix[11] = c*(u2 + v2) - w*(a*u + b*v)
                 + (w*(a*u + b*v) - c*(u2 + v2))*cosT + (a*v - b*u)*sinT;

    matrix[0] = (u2 + (v2 + w2) * cosT);
    matrix[1] = (u*v * (1 - cosT) - w*sinT);
    matrix[2] = (u*w * (1 - cosT) + v*sinT);

    matrix[4] = (u*v * (1 - cosT) + w*sinT);
    matrix[5] = (v2 + (u2 + w2) * cosT);
    matrix[6] = (v*w * (1 - cosT) - u*sinT);

    matrix[8] = (u*w * (1 - cosT) - v*sinT);
    matrix[9] = (v*w * (1 - cosT) + u*sinT);
    matrix[10] = (w2 + (u2 + v2) * cosT);

    return matrix;
}

void
destroy_rotation_matrix(double *matrix)
{
    free(matrix);
}

void
apply_rotation_matrix_to_coord_p(coord_3D *c, double *matrix)
{
    coord_3D temp = *c;

    c->i = matrix[0]*temp.i + matrix[1]*temp.j + matrix[2]*temp.k + matrix[3];
    c->j = matrix[4]*temp.i + matrix[5]*temp.j + matrix[6]*temp.k + matrix[7];
    c->k = matrix[8]*temp.i + matrix[9]*temp.j + matrix[10]*temp.k + matrix[11];
}

void 
rotate_coords_about_axis_dp_list(coord_3D **coords, int num_coords, coord_3D point, 
              vectormag_3D direction, double theta)
{
    int i;
    double *matrix;

    matrix = create_rotation_matrix(point, direction, theta);
    for(i = 0; i < num_coords; i++)
        apply_rotation_matrix_to_coord_p(coords[i], matrix);

    destroy_rotation_matrix(matrix);
}

coord_3D
get_cartesian_point_from_internal_coords(coord_3D a, coord_3D b, coord_3D c,
                      double theta, double phi, double distance)
{
    vectormag_3D lmn_x, lmn_y, lmn_z;
    double x_p, y_p, z_p;

    vectormag_3D cb = get_vector_from_coords(c, b);
    vectormag_3D ba = get_vector_from_coords(b, a);

    lmn_y = normalize_vec(get_crossprod(ba, cb));
    lmn_z = normalize_vec(cb);
    lmn_x = get_crossprod(lmn_y, lmn_z);

    x_p = distance * sin(theta) * cos(phi);
    y_p = distance * sin(theta) * sin(phi);
    z_p = distance * cos(theta);

    return (coord_3D) {lmn_x.i*x_p + lmn_y.i*y_p + lmn_z.i*z_p + c.i,
                    lmn_x.j*x_p + lmn_y.j*y_p + lmn_z.j*z_p + c.j,
                    lmn_x.k*x_p + lmn_y.k*y_p + lmn_z.k*z_p + c.k};
}

