/* File amber_prmtop_init.c begun on or around 20080424 by BLFoley
	Purpose: initialize structure info for amber prmtop file read */
//#include "mylib.h"
#include "AMBER/amber_prmtop.h"
//#include "amber_prmtop_structs.h"

void amber_prmtop_init(amber_prmtop *P){
int a=0,b=0,c=0;
char *s,*sp,**p1,**p2,*FLG,*FOR; // strings

// Set some initial values, just 'cause
P[0].NATOM  =0; // total number of atoms 
P[0].NTYPES =0; // total number of distinct atom types
P[0].NBONH  =0; // number of bonds containing hydrogen
P[0].MBONA  =0; // number of bonds not containing hydrogen
P[0].NTHETH =0; // number of angles containing hydrogen
P[0].MTHETA =0; // number of angles not containing hydrogen
P[0].NPHIH  =0; // number of dihedrals containing hydrogen
P[0].MPHIA  =0; // number of dihedrals not containing hydrogen
P[0].NHPARM =0; // currently not used
P[0].NPARM  =0; // currently not used
P[0].NEXT   =0; // number of excluded atoms
P[0].NRES   =0; // number of residues
P[0].NBONA  =0; // MBONA + number of constraint bonds
P[0].NTHETA =0; // MTHETA + number of constraint angles
P[0].NPHIA  =0; // MPHIA + number of constraint dihedrals
P[0].NUMBND =0; // number of unique bond types
P[0].NUMANG =0; // number of unique angle types
P[0].NPTRA  =0; // number of unique dihedral types
P[0].NATYP  =0; // number of atom types in parameter file, see SOLTY below
P[0].NPHB   =0; // number of distinct 10-12 hydrogen bond pair types
P[0].IFPERT =0; // set to 1 if perturbation info is to be read in
P[0].NBPER  =0; // number of bonds to be perturbed
P[0].NGPER  =0; // number of angles to be perturbed
P[0].NDPER  =0; // number of dihedrals to be perturbed
P[0].MBPER  =0; // number of bonds with atoms completely in perturbed group
P[0].MGPER  =0; // number of angles with atoms completely in perturbed group
P[0].MDPER  =0; // number of dihedrals with atoms completely in perturbed groups
P[0].IFBOX  =0; // set to 1 if standard periodic box, 2 when truncated octahedral
P[0].NMXRS  =0; // number of atoms in the largest residue
P[0].IFCAP  =0; // set to 1 if the CAP option from edit was specified
//
P[0].nS=0;
P[0].nSS=0;
P[0].nES=0;
P[0].nFF=0;

// Set integer values for strings in AMBER_PRMTOP_FLAGS 
// Set number of chars per entry from AMBER_PRMTOP_FORMATS
// To do this...
//	see how much memory we need -- 
// 	loop through AMBER_PRMTOP_FLAGS and count fields between spaces
a=0;
b=1; // will need one more field than # of spaces
while(AMBER_PRMTOP_FLAGS[a]!='\0'){
	if(AMBER_PRMTOP_FLAGS[a]==' ') b++;
	a++;}
// check the FORMATS definition, just 'cause
a=0;
c=1; // will need one more field than # of spaces 
while(AMBER_PRMTOP_FORMATS[a]!='\0'){
	if(AMBER_PRMTOP_FORMATS[a]==' ') c++;
	a++;}
if(b!=c){
	printf("b is %d ; c is %d \n",b,c);
	printf("Mismatch in AMBER FLAGS and FORMATS as defined in amber_prmtop_structs.h. Exiting.\n"); 
	exit(1);}
//	allocate memory in the amber_prmtop structure
P[0].FLAGS=(char**)calloc(b,sizeof(char*));
P[0].FORMATS=(int*)calloc(b,sizeof(int));
// Use strtok to split the flag strings and get ascii integer variants
// save each flag as the integer value of the string
FLG=strdup(AMBER_PRMTOP_FLAGS);
P[0].FLAGS[0]=strdup(strtok(FLG," "));
for(a=1;a<b;a++){P[0].FLAGS[a]=strdup(strtok(NULL," "));}
// now parse the format strings
// save each format as the number of characters allocated for each datum
FOR=strdup(AMBER_PRMTOP_FORMATS);
p1=&FOR;
s=strdup(strtok_r(FOR," ",p1)); // get a format string
p2=&s;
sp=strdup(strtok_r(s,"aAiIfFeE",p2)); // find the format identifier
sp=strdup(strtok_r(NULL,".",p2)); // information after, but before any periods
if((strchr(sp,'S')!=NULL)||(strchr(sp,'N')!=NULL)){sp[0]='0';} // EN or ES
sscanf(sp,"%d",&P[0].FORMATS[0]);
for(a=1;a<b;a++){
	s=strdup(strtok_r(NULL," ",p1)); // get a format string
	sp=strdup(strtok_r(s,"aAiIfFeE",p2)); // find the format identifier
	sp=strdup(strtok_r(NULL,".",p2)); // after, but before any periods
	if((strchr(sp,'S')!=NULL)||(strchr(sp,'N')!=NULL)){sp[0]='0';} // EN or ES
	sscanf(sp,"%d",&P[0].FORMATS[a]); 
	}
return;
}
