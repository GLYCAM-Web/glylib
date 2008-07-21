/* File amber_prmtop.h begun on 20080622 by BLFoley
 * Purpose: header file for functions that read amber prmtop 
 * files into the GLYLIB structures
 */
#include <mylib.h> 
#include <general.h> 
#include <molecules.h>
#include <declarations.h>
#include "amber_prmtop_structs.h"

void amber_prmtop_init(amber_prmtop *P);
void read_amber_prmtop_asis(fileset F,amber_prmtop *P);
assembly parse_amber_prmtop(amber_prmtop *P);
void find_molecules_molbond_array(int nMB, molbond *MB, int NATOM, molindex *AT);

