/* File find_molecules.c begun on 20080625 by BLFoley 
	Purpose: Contain functions designed to determine which atoms and/or
	residues are parts of molecules and to separate them into molecules.

	For example, find_molecules_molbond_array will:
		Given an array of molbond structures, use the bonding to 
		determine which atoms (and residues) are parts of molecules.
		The molecules found will also be numbered. */

// change this before adding to library -- doesn't need whole prmtop header
// this is here because the first function written was written for a read
// of an amber prmtop file.
#include "AMBER/amber_prmtop.h"

typedef struct {
	int n; // number of
	int *i; // n of these 
} tempintset;

void follow_find_molecule_amber_prmtop_bond_target(molbond *MB, molindex *AT, tempintset *SI, tempintset *TI, int i, int mN);

/***************  find_molecules_molbond_array() ******************/
/* Finds molecules in a molbond array based on bonding patterns 
	Updates an array of NATOM molindices with molecule information.
	Also writes molecule info to the molbond array.  
	In both cases, it overwrites existing information. 
	NOTES:  
	1.  There should be no useful molecule or residue information
	in the molbond array (it will be igored).
	2.  The index, i, in the molindex structures should correspond to the
	absolute atom number as indexed in the molindex array.
	This way, any residue and atom info you have previously will
	be retained -- and it makes this function loads easier to write.
	Later, you can renumber your residues and atoms however you like,
	but any existing molecule information will be overwritten.
*/
void find_molecules_molbond_array(int nMB, molbond *MB, int NATOM, molindex *AT){
int a=0,currmol=0,localmol=0; // current and local molecule numbers
tempintset *SI,*TI;

// set all AT molecule indices to -1 (not seen)
for(a=0;a<NATOM;a++){
//printf("(find mols) AT[%d].a is %d ; .i is %d \n",a,AT[a].a,AT[a].i);
	AT[a].m=-1;
	}

// index each source and target to all relevant indices in the AT array
//	calloc space  (other way???)
SI=(tempintset*)calloc(NATOM,sizeof(tempintset));
TI=(tempintset*)calloc(NATOM,sizeof(tempintset));
for(a=0;a<NATOM;a++){ // allocate space 
	SI[a].i=(int*)calloc(1,sizeof(int));
	TI[a].i=(int*)calloc(1,sizeof(int)); }
// set backward indices for the atoms to the bonds
//for(a=0;a<nMB;a++){
//printf("MB[%d] : s = %d ; t = %d\n",a,MB[a].s.i,MB[a].t.i);}
for(a=0;a<nMB;a++){ 
//printf("SI[MB[a].s.i].n is %d TI[MB[a].t.i].n is %d \n",SI[MB[a].s.i].n,TI[MB[a].t.i].n);
	SI[MB[a].s.i].n++;
	TI[MB[a].t.i].n++;
	SI[MB[a].s.i].i=(int*)realloc(SI[MB[a].s.i].i,SI[MB[a].s.i].n*sizeof(int));
	TI[MB[a].t.i].i=(int*)realloc(TI[MB[a].t.i].i,TI[MB[a].t.i].n*sizeof(int));
	SI[MB[a].s.i].i[SI[MB[a].s.i].n-1]=a;
	TI[MB[a].t.i].i[TI[MB[a].t.i].n-1]=a;
//printf("  MB[%d].s.i is %d ; SI[MB[a].s.i].n-1 is %d \n",a,MB[a].s.i,SI[MB[a].s.i].n-1);
//printf("  MB[%d].t.i is %d ; TI[MB[a].t.i].n-1 is %d \n",a,MB[a].t.i,TI[MB[a].t.i].n-1);
	}
//for(a=0;a<NATOM;a++){
//printf("SI[%d].n is %d\n",a,SI[a].n); 
//for(b=0;b<SI[a].n;b++){
//printf("\tSI[%d].i[%d] is %d\n",a,b,SI[a].i[b]);
//}
//printf("TI[%d].n is %d\n",a,TI[a].n); 
//for(b=0;b<TI[a].n;b++){
//printf("\tTI[%d].i[%d] is %d\n",a,b,TI[a].i[b]);
//}
//}
//for(a=0;a<nMB;a++){
//printf("MB[%d] : s = %d ; t = %d\n",a,MB[a].s.i,MB[a].t.i);}

// Find the molecules
for(a=0;a<nMB;a++){ // loop through the molbond array 
	localmol=-1; // = we don't know the local molecule number yet
	// see if s or t already belong to a molecule 
	// if so, set localmol to that number
//printf("MB[%d].s.i is %d ;  AT[MB[a].s.i].m is %d  ; MB[a].t.i is %d ;  AT[MB[a].t.i].m is %d  \n",a,MB[a].s.i,AT[MB[a].s.i].m,MB[a].t.i,AT[MB[a].t.i].m);
	if(AT[MB[a].s.i].m!=-1) {localmol=AT[MB[a].s.i].m;}
	if(AT[MB[a].t.i].m!=-1) {
		if(localmol!=-1){
			if(AT[MB[a].t.i].m!=AT[MB[a].s.i].m){
			mywhine("molecule mismatch in find_molecules_molbond_array");} }
		if(localmol==-1){ localmol=AT[MB[a].s.i].m=AT[MB[a].t.i].m;}
		else {AT[MB[a].t.i].m=localmol;}
		}
	// if we still don't know the local molecule number
	if(localmol==-1){ // set localmol to currmol and increment currmol
		localmol=currmol;
		currmol++;
		// set s and t AT's as belonging to localmol
		AT[MB[a].s.i].m=AT[MB[a].t.i].m=localmol;
		}
	// set molecule info local to the bonds, too
	MB[a].s.m=MB[a].t.m=localmol; 
	follow_find_molecule_amber_prmtop_bond_target(MB,AT,SI,TI,MB[a].t.i,localmol);
	}
// check our work:
for(a=0;a<nMB;a++){ 
	if(MB[a].s.m==-1) {mywhine("source MB molecule not defined at end of assignment.");}
	if(MB[a].t.m==-1) {mywhine("target MB molecule not defined at end of assignment.");}
	}
for(a=0;a>NATOM;a++){ // assign lone atoms to molecules
	if(AT[a].m==-1){
		AT[a].m=currmol;
		currmol++;
		} 
	}

return;
}

/************** follow_find_molecule_amber_prmtop_bond_target() ************/
/* Follow bonds along the molecule recursively */
void follow_find_molecule_amber_prmtop_bond_target(molbond *MB, molindex *AT, tempintset *SI, tempintset *TI, int i, int mN){
int b=0; // for counting 

//Each bond indicator has a source atom and a target atom.  So, for each target atom,
//(1) find all bonds for which that target is a source
//printf("in recursive follow at top\n");
//printf("m is %d ; the atom index, i, is %d ; SI[i].n is %d ; TI[i].n is %d \n",mN,i,SI[i].n,TI[i].n);

for(b=0;b<SI[i].n;b++){
//printf("SI-i=%d (recursive follow, m=%d) the first source bond is %d  \n",i,mN,SI[i].i[b]);
//printf("\t the source-target atoms are %d-%d ; its molecule is %d \n",MB[SI[i].i[b]].s.i,MB[SI[i].i[b]].t.i,AT[MB[SI[i].i[b]].s.i].m);
	if((AT[MB[SI[i].i[b]].s.i].m!=-1)&&(AT[MB[SI[i].i[b]].s.i].m!=mN)){mywhine("AT[MB[x].s.i].m!=-1)&&(AT[MB[x].s.i].m!=mN, in recursive follow\n");}
	//if(AT[MB[SI[i].i[b]].t.i].m!=-1) continue;  
	AT[MB[SI[i].i[b]].s.i].m=mN; // set molecule for the source
	if(AT[MB[SI[i].i[b]].t.i].m==-1){
		AT[MB[SI[i].i[b]].t.i].m=mN; // set molecule for the target
		follow_find_molecule_amber_prmtop_bond_target(MB,AT,SI,TI,MB[SI[i].i[b]].t.i,mN); }
	}
	//-- set their targets as belonging to the molecule
		//(check to be sure not already marked someone else's)
	//-- call this function for each of the targets
//(2) find all bonds for which that target is a target
	//-- set their sources as belonging to the molecule
	//-- call this function for each of the sources
// best to do this recursively....
// also set all targets for that atom as belonging
for(b=0;b<TI[i].n;b++){
//printf("TI-i=%d (recursive follow, m=%d) the first source bond is %d  \n",i,mN,TI[i].i[b]);
//printf("\t the source-target atoms are %d-%d ; its molecule is %d \n",MB[TI[i].i[b]].t.i,MB[TI[i].i[b]].s.i,AT[MB[TI[i].i[b]].t.i].m);
	if((AT[MB[TI[i].i[b]].t.i].m!=-1)&&(AT[MB[TI[i].i[b]].t.i].m!=mN)){mywhine("AT[MB[x].t.i].m!=-1)&&(AT[MB[x].t.i].m!=mN, in recursive follow\n");}
	//if(AT[MB[TI[i].i[b]].t.i].m!=-1) continue;
	AT[MB[TI[i].i[b]].t.i].m=mN; // set molecule for the target
	if(AT[MB[TI[i].i[b]].s.i].m==-1){ // set molecule for the target
		AT[MB[TI[i].i[b]].s.i].m=mN; // set molecule for the source
		follow_find_molecule_amber_prmtop_bond_target(MB,AT,SI,TI,MB[TI[i].i[b]].s.i,mN); }
	}
return;
}
