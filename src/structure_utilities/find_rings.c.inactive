/** \file find_rings.c
 *  \brief Functions for finding rings.
 *
 *  begin 20110421 BLFoley 
 */

/**  set_smallest_rings_from_residue_atom_nodes finds rings in a residue.
 *
 *  This function might not work properly unless the bonding trees were 
 *  set according to the internal "from bonds" functions.  However, even 
 *  if the function thinks this might be the case, it tries anyway.
 *
 *  If there is position information in either the main coordinate or 
 *  the first alternate coordinate set, the function will set the total
 *  length (s) based on that information.  If position information is 
 *  not available, it will be set to zero.
 */
#include <stdlib.h>
#include <molecules.h>

void set_smallest_rings_from_residue_atom_nodes(residue *r)
{
int nri=0,ni=0,no=0,ai=0,bi=0,si=0,in=0,nrings=0;
int itmp=0,nr_localmax=0,*start,*nstarts,nrstarts=0,niimax=0;
int ni1=0,ni2=0;
char isorigin='n';
molindex_set *RING;
struct
    {
    int *steps;
    ensindex *parent;
    } path_follower;
path_follower *pf;

struct
	{
	int *steps;
	int *child;
	}
ring_follower *rf;

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

isorigin='n';
for(ai=0;ai<r[0].na;ai++)
    {
    if(r[0].aT[ai].ni<0){mywhine("Unexpected value of r[0].aT[ai].ni in set_smallest_rings_from_residue_atom_nodes.");}
    if(r[0].aT[ai].no<0){mywhine("Unexpected value of r[0].aT[ai].no in set_smallest_rings_from_residue_atom_nodes.");}
    if(r[0].aT[ai].isorigin=='Y') isorigin = 'y';
    }
if(isorigin!='y')
	{
	printf("\nWARNING: No atom in tree set as origin in \n\tset_smallest_rings_from_residue_atom_nodes.\n");
	printf("This will almost certainly cause trouble.\n");
	printf("Setting first atom as origin and trying anyway.\n");
	printf("If other things fall apart, this is a possible reason.\n");
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
        Record lots of other bits of information.
*/
nrings=0;
nstarts=(int*)calloc(r[0].na,sizeof(int));
for(ai=0;ai<r[0].na;ai++) 
    { 
    if(r[0].aT[ai].isorigin=='Y') nstarts[0]= ( r[0].aT[0].ni+1 ) * ( r[0].aT[0].ni ) / 2;
    else nstarts[ai] = ( r[0].aT[ai].ni-1 ) * ( r[0].aT[ai].ni-2 ) / 2; 
    if(nstarts[ai]>0)
        {
        nrstarts++;
        nrings += nstarts[ai];
        if(nstarts[ai]>nr_localmax) nr_localmax=nstarts[ai];
        if(r[0].aT[ai].ni>niimax) niimax=r[0].aT[ai].ni;
        }
    }
start=(int*)calloc(nrstarts,sizeof(int));
nri=0;
si=0;
for(ai=0;ai<r[0].na;ai++)
    { 
    if(nstarts[ai]>0)
        {
	start[si]=ai;
	si++;
	}
    nri+=nstarts[ai]; /* a quick code check */ 
    }
if(nri!=nrings){mywhine("nri!=nrings in set_smallest_rings_from_residue_atom_nodes.");}

/*
        Assign accounting spaces.

*/
rf = (ring_follower*) calloc(r[0].na, sizeof(ring_follower));
for(ai=0;ai<r[0].na;ai++)
    { 
    rf[ai].steps=(int*)calloc(niimax,sizeof(int));
    rf[ai].child=(int*)calloc(niimax,sizeof(int));
    }
RING = (molindex_set*) calloc(nrings, sizeof(molindex_set));

/*
    1.  For each atom with multiple incoming bonds.
*/

oh.... change this to be recursive.  Much simpler.  
	for each incoming, save only shortest path to origin




for(si=0;si<nrstarts;si++)
    {
/*  
    1.0  Initialize the accounting spaces.
*/ 
    for(ai=0;ai<r[0].na;ai++) {
    for(nii=0;nii<niimax;nii++) {
        rf[ai].steps[nii] = rf[ai].child[nii] = -1;
        } } 
/*
    1.1  Follow each incoming up to the origin.
        Label them all in temporary rings.
*/ 
    for(nii=0;nii<r[0].aT[starts[si]].ni;nii++)
    	{
	steps=0;
        isorigin='n';
	nri=si;
        while(isorigin=='n')
		{
		rf[nri].steps[nii] = (steps++);
		rf[nri].child[nii] = r[0].aT[nri].i[0][some variable...];
		}
	}
/*
    1.2  Find where the two incoming sets intersect.

        If more than one intersection, find smallest.
        If two or more are the same size, increase nrings.
	Mark the intersection (or multiple intersections)
        that should be followed.
*/
    ni1=0;
    ni2=1;
    for(nri=0;nri<nstarts[ai];nri++)
        {
        ni2++;
        if(ni2==r[0].aT[starts[si]].ni)
            {
            ni1++;
            ni2=ni1+1;
            }
        }
/*
    1.3  Add the nodes to a "ring".
*/
    }

return;
}
