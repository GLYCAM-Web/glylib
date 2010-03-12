/** \file amber_prmtop.h Header file for functions that read amber prmtop 
 * files into the GLYLIB structures

begun on 20080622 by BLFoley
 */

#if !defined(GLYLIB_AMBER_PRMTOP_HEADER)
#define GLYLIB_AMBER_PRMTOP_HEADER

#include <mylib.h> 
#include <general.h> 
#include <molecules.h>
#include <declarations.h>
#include <AMBER/amber_prmtop_structs.h>

void amber_prmtop_init(amber_prmtop *P);
void read_amber_prmtop_asis(fileset F,amber_prmtop *P);
assembly parse_amber_prmtop(amber_prmtop *P);
void find_molecules_molbond_array(int nMB, molbond *MB, int NATOM, molindex *AT);
assembly load_amber_prmtop(fileset F);
void add_trajcrds_to_prmtop_assembly(
        fileset F, ///< The trajectory file containing the data *Must* be open
        assembly *A, ///< Assembly ("incoming") to which the crds should be added
        char ftype, ///< 'c' if coordinate ; 'v' if velocity ; 'r' if restart (might have both)
        int offset ///< read in the offset-th trajectory -- starts with zero=current
        );
void deallocateAmberPrmtopSection(amber_prmtop_section *aps);
void deallocateAmberPrmtop(amber_prmtop *ap);

#endif
