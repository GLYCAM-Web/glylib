/** file PDB.h 

    Contains header info relevant to reading and writing PDB files. */

#if !defined(GLYLIB_PDB)
#define GLYLIB_PDB
#include <molecules.h>
#include <gly_codeutils.h>
#include <gly_fileutils.h>
#include <structures.h>

fileslurp get_ensemble_PDB_ATOM_lines(ensemble *E,char isource, int savei,char raltname, int xs);

fileslurp get_assembly_PDB_ATOM_lines(assembly *A,char isource, int savei,char raltname, int xs);

const char *get_PDB_line_for_ATOM(atom *a, residue *r,int ai, int ri, int asave, char raltname, int xs);

void make_ATOM_HETATM(char *pdbline);

fileslurp get_residue_PDB_ATOM_lines(residue *r,int ri,int ainit, int rsave, int asave,char raltname, int xs);

fileslurp get_molecule_PDB_ATOM_lines(molecule *mol,int rinit,int ainit, int rsave,int asave, char oneres, char raltname, int xs);

void outputMolPDB(molecule*,char*);  /* Writes a pdb using a given molecule -- deprecated */
void outputAsmblPDB(assembly*,char*);/* Writes a pdb using a given assembly -- deprecated */

#endif
