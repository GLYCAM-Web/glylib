#include "molecules.h"
#include <stdlib.h>

coord **
atoms_to_coord_list(atom **atoms, int num_atoms)
{
    int i;

    coord **coords = (coord **)malloc(num_atoms * sizeof(coord *));
    for(i=0; i<num_atoms; i++)
        coords[i] = &atoms[i]->x;

    return coords;
}

