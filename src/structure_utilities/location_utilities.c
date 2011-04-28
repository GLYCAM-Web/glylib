/** \file location_utilities.c
*   Purpose:  Functions for finding atoms and molecules
*   Begin 20110428 BLFoley
*
*/
#include <mylib.h>
#include <molecules.h>

/*
Find atoms named N -- searches for a.N
*/
molindex_set find_residue_atoms_by_N(residue *r, const char *name)
{
  molindex_set moli_set;
  int i,match*,matches=0;

  match=(int*)calloc(r[0].na,sizeof (int)); 
  for(i=0;i<r[0].na;i++)
    {
       if(strcmp(r[0].N,name)==0)
           {
           match[i]=1;
           matches++;
           }
    } 
  moli_set.nP=matches;
  moli_set.P=(molindex*)calloc(matches,sizeof(molindex));
  matches=0;
  for(i=0;i<r[0].na;i++)
    {
    if(match[i]==1)
       {
       moli_set.P[matches]=r[0].a[i].moli;
       matches++;
       }
    if(matches==moli_set.nP){break;}
    }
return moli_set;
}

molindex_set find_molecule_atoms_by_N(molecule *m, const char *name)
{
}
molindex_set find_molecule_residues_by_N(molecule *m, const char *name)
{
}
/*
Find atoms numbered n -- searches for a.n
*/
molindex_set find_residue_atoms_by_n(residue *r, int number)
{
}
molindex_set find_molecule_atoms_by_n(molecule *m, int number)
{
}
molindex_set find_molecule_residues_by_n(molecule *m, int number)
{
}
