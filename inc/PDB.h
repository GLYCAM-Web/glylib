/** file PDB.h 

    Contains header info relevant to reading and writing PDB files. */

#if !defined(GLYLIB_PDB)
#define GLYLIB_PDB
#include <molecules.h>
#include <gly_codeutils.h>
#include <gly_fileutils.h>
#include <structures.h>

/*
    PDB Writing Utilities

    In the following:

        isource = 'n' to use values stored in r.n and a.n
                  'i' to assign numbers automatically
        ai = the index to use for a current atom (serial)
        ri = the index to use for a current residue (resSeq)
        ainit = the atom number (serial) to start a list with
                if -1 then the value saved in a.n will be used
        rinit = the residue number (resSeq) to start a list with
                if -1 then the value saved in r.n will be used
        isave = common setting for asave and rsave
        asave = the index in a.i where the assigned serial should be saved
                this is used for setting CONECT and LINK cards 
                if -1, will not be saved
        rsave = the index in r.i where the assigned resSeq should be saved
                this is used for setting CONECT and LINK cards 
                if -1, will not be saved 
        raltname = if 'y', use the residue name stored in r.altname
                   instead of r.N -- if used with oneres='y', set them
                   all to be the same
        oneres = 'y' to make the whole molecule one residue
                 'n' to leave it as separate residues
        
*/
int NUMAT;  /* Would like to get rid of this eventually */

/** \addtogroup PDB
 * @{
 */
fileslurp get_ensemble_PDB_ATOM_lines(ensemble *E,char isource, int savei,char raltname, int xs);

fileslurp get_assembly_PDB_ATOM_lines(assembly *A,char isource, int savei,char raltname, int xs);
fileslurp get_assembly_PDB_CONECT_lines(assembly *A, char isource, int savei);
fileslurp get_molecule_PDB_CONECT_lines(molecule *m, int savei);
get_atom_PDB_CONECT_lines_assembly(assembly *A, atom *a, int savei);

const char *get_PDB_line_for_ATOM(atom *a, residue *r,int ai, int ri, int asave, char raltname, int xs);

void make_ATOM_HETATM(char *pdbline);

fileslurp get_residue_PDB_ATOM_lines(residue *r,int ri,int ainit, int rsave, int asave,char raltname, int xs);

fileslurp get_molecule_PDB_ATOM_lines(molecule *mol,int rinit,int ainit, int rsave,int asave, char oneres, char raltname, int xs);

void outputMolPDB(molecule*,char*);  /* Writes a pdb using a given molecule -- deprecated */
void outputAsmblPDB(assembly*,char*);/* Writes a pdb using a given assembly -- deprecated */
/** @}*/

#endif
