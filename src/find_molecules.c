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
int a=0,b=0,currmol=0,localmol=0; // current and local molecule numbers
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
//printf("  MB[%d].s.i is %d ; SI[MB[a].s.i].n-1 is %d \t",a,MB[a].s.i,SI[MB[a].s.i].n-1);
//printf("  MB[%d].t.i is %d ; TI[MB[a].t.i].n-1 is %d \n",a,MB[a].t.i,TI[MB[a].t.i].n-1);
	}
for(a=0;a<NATOM;a++){
//printf("SI[%d].n is %d\n",a,SI[a].n); 
for(b=0;b<SI[a].n;b++){
//printf("\tSI[%d].i[%d] is %d\n",a,b,SI[a].i[b]);
}
//printf("TI[%d].n is %d\n",a,TI[a].n); 
for(b=0;b<TI[a].n;b++){
//printf("\tTI[%d].i[%d] is %d\n",a,b,TI[a].i[b]);
}
}
for(a=0;a<nMB;a++){
//printf("MB[%d] : s = %d ; t = %d\n",a,MB[a].s.i,MB[a].t.i);
}

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
			mywhine("molecule mismatch in find_molecules_molbond_array");} 
			}
		if(localmol==-1){ localmol=AT[MB[a].s.i].m=AT[MB[a].t.i].m;}
		else {AT[MB[a].t.i].m=localmol;}
		}
	// if we still don't know the local molecule number
	if(localmol==-1){ // set localmol to currmol and increment currmol
		localmol=currmol;
		currmol++;
		// set s and t AT's as belonging to localmol
		//AT[MB[a].s.i].m=AT[MB[a].t.i].m=localmol;
		}
	// set molecule info local to the bonds, too
	//MB[a].s.m=MB[a].t.m=localmol; 
	//follow_find_molecule_amber_prmtop_bond_target(MB,AT,SI,TI,MB[a].t.i,localmol);
	// try this way instead... Start with sources only
	follow_find_molecule_amber_prmtop_bond_target(MB,AT,SI,TI,MB[a].s.i,localmol);
	}
// check our work:
for(a=0;a<nMB;a++){ 
//printf("a=%d ; MB[a].s.m=%d ; MB[a].t.m=%d\n",a,MB[a].s.m,MB[a].t.m);
	if(MB[a].s.m==-1) {mywhine("source MB molecule not defined at end of assignment.");}
	if(MB[a].t.m==-1) {mywhine("target MB molecule not defined at end of assignment.");}
	}
for(a=0;a>NATOM;a++){ // assign lone atoms to molecules
	if(AT[a].m==-1){
		AT[a].m=currmol;
		currmol++;
		} 
	}

for(a=0;a<NATOM;a++){ // allocate space 
	free(SI[a].i);
	free(TI[a].i);
	}
free(SI);
free(TI);
return;
}

/************** follow_find_molecule_amber_prmtop_bond_target() ************/
/* Follow bonds along the molecule recursively */
void follow_find_molecule_amber_prmtop_bond_target(molbond *MB, molindex *AT, tempintset *SI, tempintset *TI, int i, int mN){
int b=0,c=0; // for counting 
// for sanity:
int 	sourcei=0,  // index for the instance of this atom being a bond source
	targeti=0,  // index for the instance of this atom being a bond target
	sbondi=0,   // index into the molbond array for this source index
	tbondi=0,   // index into the molbond array for this target index
	stbondi=0;  // index into molbond array for sources of targets or targets as sources

//Each bond indicator has a source atom and a target atom.  So, for each target atom,
//(1) find all bonds for which that target is a source
//printf("in recursive follow at top\n");
//printf("m is %d ; the atom index, i, is %d ; SI[i].n is %d ; TI[i].n is %d \n",mN,i,SI[i].n,TI[i].n);

for(b=0;b<SI[i].n;b++){
//printf("SI-i=%d (recursive follow, m=%d) the first source bond is %d  \n",i,mN,SI[i].i[b]);
//printf("\t the source-target atoms are %d-%d ; its molecule is %d \n",MB[SI[i].i[b]].s.i,MB[SI[i].i[b]].t.i,AT[MB[SI[i].i[b]].s.i].m);
	sourcei=SI[i].i[b]; // the b-th target associated with this source
	sbondi=MB[sourcei].s.i; // the bonding index info associated with this source
	if((AT[sbondi].m!=-1)&&(AT[sbondi].m!=mN)){mywhine("AT[sbondi].m!=-1)&&(AT[sbondi].m!=mN, in SI recursive follow\n");}
	AT[sbondi].m=MB[sourcei].s.m=MB[sourcei].t.m=mN; // set molecule for the source
//printf("sbondi=%d ; just set MB[%d].s/t.m=%d\n",sbondi,sourcei,MB[sourcei].s.m);
	tbondi=MB[sourcei].t.i;
	if((AT[tbondi].m!=-1)&&(AT[tbondi].m!=mN)){mywhine("AT[tbondi].m!=-1)&&(AT[tbondi].m!=mN, in SI recursive follow\n");} 
	if(AT[tbondi].m==-1){
		AT[tbondi].m=MB[sourcei].t.m=mN; // set molecule for the target
		for(c=0;c<SI[tbondi].n;c++){ // for each other bond for which this target is a source
			targeti=SI[tbondi].i[c];
			stbondi=MB[targeti].s.i;
			// follow each of the sources in those other targets, but only if not already assigned
			if((AT[stbondi].m!=-1)&&(AT[stbondi].m!=mN)){
				mywhine("AT[stbondi].m!=-1)&&(AT[stbondi].m!=mN, in SI-TI recursive follow\n");}
			follow_find_molecule_amber_prmtop_bond_target(MB,AT,SI,TI,MB[targeti].s.i,mN); 
			}
//printf("just set MB[%d].t.m=%d\n",sourcei,MB[sourcei].t.m);
		for(c=0;c<TI[tbondi].n;c++){ // for each other bond for which this target is a target
			targeti=TI[tbondi].i[c];
			stbondi=MB[targeti].s.i;
			// follow each of the sources in those other targets, but only if not already assigned
			if((AT[stbondi].m!=-1)&&(AT[stbondi].m!=mN)){
				mywhine("AT[stbondi].m!=-1)&&(AT[stbondi].m!=mN, in SI-TI recursive follow\n");}
			follow_find_molecule_amber_prmtop_bond_target(MB,AT,SI,TI,MB[targeti].s.i,mN); 
			}
		}
	}
// If this source is also a target for another bond, follow that source, too
for(b=0;b<TI[i].n;b++){
//printf("TI-i=%d (recursive follow, m=%d) the first source bond is %d  \n",i,mN,TI[i].i[b]);
//printf("\t the target-source atoms are %d-%d ; its molecule is %d \n",MB[TI[i].i[b]].t.i,MB[TI[i].i[b]].s.i,AT[MB[TI[i].i[b]].t.i].m);
	targeti=TI[i].i[b];
	tbondi=MB[targeti].t.i;
	if(AT[tbondi].m!=mN){mywhine("AT[tbondi].m!=mN in TI recursive follow\n");} // should already be set...
	//if((AT[tbondi].m!=-1)&&(AT[tbondi].m!=mN)){mywhine("AT[tbondi].m!=-1)&&(AT[tbondi].m!=mN, in TI recursive follow\n");}
	//AT[MB[TI[i].i[b]].t.i].m=mN; // set molecule for the target
	sbondi=MB[TI[i].i[b]].s.i;
	if(AT[sbondi].m==-1){ // if this source has not already been seen
		//AT[sbondi].m=mN; // set molecule for the source -- should not be necessary
		follow_find_molecule_amber_prmtop_bond_target(MB,AT,SI,TI,sbondi,mN); }
	}

return;
}
