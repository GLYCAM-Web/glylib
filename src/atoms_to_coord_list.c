#include <stdlib.h>
#include <molecules.h>

coord_3D **
atoms_to_coord_list(atom **atoms, int num_atoms)
{
    int i;

    coord_3D **coords = (coord_3D **)malloc(num_atoms * sizeof(coord_3D *));
    for(i=0; i<num_atoms; i++)
        coords[i] = &atoms[i]->x;

    return coords;
}

