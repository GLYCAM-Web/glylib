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
void set_rings_from_residue_atom_nodes(residue *r)
{

/*
    0.  Check that the nodes are set and seem sane.
        For now, check that there is only one ring.
        Add multi-ring capability later.
        If the origin has one or more incoming, there 
        might be a problem.
*/
/*
    1.  Find multiple incoming bonds at one atom.
*/
/*
    2.  Follow both incoming up to the origin.
*/
/*
    3.  Find where the two incoming sets intersect.

        If more than one intersection, find smallest.
        If more than one intersection, why did it pass
        test 0?  (until multi-ring is written)
*/
/*
    4.  
*/
return;
}
