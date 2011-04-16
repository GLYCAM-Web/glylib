/* Author: Michael Tessier
	Extended and updated by BLFoley */
#include <mylib.h>
#include <PDB.h>
#include <gly_codeutils.h>

fileslurp get_molecule_PDB_CONECT_lines(molecule *m,int rsave,int asave);

fileslurp 
get_molecule_PDB_CONECT_lines(molecule *m,int rsave,int asave)
{
fileslurp FM;


/* 
    A couple brief sanity checks.
*/
if((m[0].r[0].ni<=rsave)||(m[0].r[0].a[0].ni<=asave))
    {
    mywhine("Atom serial location specified does not seem to exist in get_molecule_PDB_CONECT_lines.");
    }

printf("get_molecule_PDB_CONECT_lines isn't complete yet...  please complete it!");
exit(0);

return FM;
}
