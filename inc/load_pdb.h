#include <molecules.h>
#include <structures.h>
#include <mylib.h>
#define LA (int)0

assembly* load_pdb(char* file_name);
int howManyMolecules();
assembly* getAssembly();
int findTotalResidue(int start);
int endOfMol(linedef* line);
int isAtom(linedef* line);
void getResInfo(residue* res, int start);
void init_struct();
linetype get_type(char*);
void pdb_def(); /* contains definitions of lines and fields */
// moving this to mylib where the other read error functions are
//void read_neek(const char*,int,int);
void rwm_line(int); /* read in pdb, call to modify, write out pdb */

FILE *IN, *OUT, *OUTC;
linedef line[1];
int WARN,SILENT; /* default is "false" (=1). "true" is 0 */
char froot[200],suf[6],sufc[20];
char ACT; // Flag tells program what to do
int INWC,DATE; /* number of lines in input file, date */
int UNCUT,UNx,UNy,UNz; /* flags for uncutting the box */
int LASTRES,LASTOKX,LASTOKY,LASTOKZ;
float UNCTOL,CRYX,CRYY,CRYZ,LASTX,LASTY,LASTZ;
int DEBUG;
