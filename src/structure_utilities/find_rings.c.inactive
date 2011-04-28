/** \file find_rings.c
 *  \brief Functions for finding rings.
 *
 *  begin 20110421 BLFoley 
 */

/*  The following function should soon be improved.  At the moment, it is 
 *  only being written for a single ring in a single residue.  The author
 *  pleads overextension.
 *
 *  This function will not work properly unless the bonding trees were 
 *  set according to the internal "from bonds" functions.
 */
#include <stdlib.h>
#include <molecules.h>

void set_smallest_rings_from_residue_atom_nodes(residue *r)
{
int nri=0,ni=0,no=0,ai=0,bi=0,si=0,nrings=0;
int *label[2],*steps[2],*child[2],*start;
struct
    {
    int np;
    molindex *p;
    ringset *next;
    } ringset;
ringset rings, *rcurrent;

/*
    0.  Check that the nodes are set and seem sane.
        If the origin has one or more incoming, there 
        might be a problem.
*/

if(r[0].aT==NULL){mywhine("r[0].aT==NULL in set_smallest_rings_from_residue_atom_nodes.");}
if(r[0].aT[0].isorigin!='Y'){mywhine("r[0].aT[0] is not origin in set_smallest_rings_from_residue_atom_nodes.");}

if(r[0].na<3)
    {  /* If there aren't enough atoms to make a ring. */
    r[0].nring=0;
    return;
    }

for(ai=0;ai<r[0].na;ai++)
    {
    if(r[0].aT[ai].ni<0){mywhine("Unexpected value of r[0].aT[ai].ni in set_smallest_rings_from_residue_atom_nodes.");}
    if(r[0].aT[ai].no<0){mywhine("Unexpected value of r[0].aT[ai].no in set_smallest_rings_from_residue_atom_nodes.");}
    }

/*
        If there are incoming bonds to the first atom, the atom nodes probably 
        aren't set properly for this function.  But, give it a go anyhow.  
*/
if(r[0].aT[0].ni>0)
    {
    printf("\nWARNING: r[0].aT[0].ni>0 in set_smallest_rings_from_residue_atom_nodes.\n");
    printf("So,the atom_node bonding probably isn't right for this function.\n");
    printf("Ignoring that and forging on anyway.\n");
    }
if(r[0].aT[0].isorigin!='Y')
    {
    printf("\nWARNING: r[0].aT[0] is not origin in set_smallest_rings_from_residue_atom_nodes.\n");
    printf("So,the atom_node bonding probably isn't right for this function.\n");
    printf("Ignoring that and forging on anyway.\n");
    }

/*
        Check number of rings.  Record start locations.
*/
nrings=0;
nrings+=r[0].at[0].ni;
for(ai=1;ai<r[0].na;ai++) { nrings += ( r[0].at[ai].ni-1 ); }
start=(int*)calloc(nrings,sizeof(int));
nri=0;
if(r[0].at[0].ni>0)
    {
    start[0]=0;
    nri++;
    }
for(ai=1;ai<r[0].na;ai++)
    {
    for(bi=0;bi<r[0].at[ai].ni-1;bi++)
        {
        start[nri]=ai;
        nri++;
        }
    }
if(nri!=nrings){mywhine("nri!=nrings in set_smallest_rings_from_residue_atom_nodes.");}

/*
        Assign accounting spaces.
*/
label[0] = (int*)malloc(r[0].na,sizeof(int));
steps[0] = (int*)malloc(r[0].na,sizeof(int));
child[0] = (int*)malloc(r[0].na,sizeof(int));
label[1] = (int*)malloc(r[0].na,sizeof(int));
steps[1] = (int*)malloc(r[0].na,sizeof(int));
child[1] = (int*)malloc(r[0].na,sizeof(int));

/*
    1.  For each atom with multiple incoming bonds.
*/
for(si=0;si<nri;si++)
    {
/*  
    1.0  Initialize the accounting spaces.
*/
    for(bi=0;bi<2;bi++)
        {
        for(ai=0;ai<r[0].na;ai++)
            {
            label[bi][ai]=steps[bi][ai]=child[bi][ai]=-1;
            }
        }

/*
    1.1  Follow each incoming up to the origin.
        Label them all in temporary rings.
*/
    

/*
    1.2  Find where the two incoming sets intersect.

        If more than one intersection, find smallest.
        If two or more are the same size, increase nrings.
	Mark the intersection (or multiple intersections)
        that should be followed.
*/
/*
    1.3  Add the nodes to a "ring".
*/
    }

return;
}
