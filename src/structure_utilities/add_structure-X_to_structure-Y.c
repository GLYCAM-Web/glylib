#include <mylib.h>
#include <molecules.h>


/** \file add_structure-X_to_structure-Y.c Functions for combining structures.

Begun on 20100203 by BLFoley and expanded, as per usual, on a crisis to 
crisis basis...
*/

/******************** add_assembly_to_ensemble() *********************/
/* Function add_assembly_to_ensemble(assembly *A, ensemble *E) will add
the contents of the assembly pointed to by A to the ensemble pointed to
by E.  Assuming the internal accounting in E has been done properly, the
contents of A will not overwrite any preexisting data in E.

The function will set assembly pointers from within the ensemble.  

Please carefully read the notes below.

This function does NOT update molecular weight, COM, etc., because the 
programmer needs to specify which mass info, etc., should be used.  You 
have to use other functions to set those values.

*** NOTE ***  A should *NOT* already be a member of E.  If A is already
a member of E, unexpected things might happen in your programs.  See the
next note for more info.

*** NOTE ***  If you use some sort of weird numbering in your molindex
elements, you need to copy those elsewhere (say to an ensi).  This function
will (re)set the moli's to be what it needs.

*** NOTE ***  The **a and **r entries in A should ONLY BE POINTING to data
in the normal m.r.a structures in **m.  Any data that has copied to 
those locations might not be freed, which is bad for your data management.
(risky to free here without knowing -- can't free empty space.)

*** A VERY IMPORTANT NOTE ***  

	DO NOT FREE THE ASSEMBLY.  You can reset the assembly pointer to 
	go somewhere else, but do not free the memory.  If you do, you will
	be freeing part of the ensemble.
	
	Here is what the function does:

	1. Allocates and/or reallocates E so that it can accomodate the
		information found in A.
	2. Copies atom information in A into the relevant locations in E.
		* If A's information is not arranged in M.r.a format, the
		function will currently complain and exit.  
	3. Make a list of the molecule indices for each atom in the double-
		pointed atom and residue lists of A.  Also records the offset
		between A's m.r.a sets and E's.  The molecules are all added
		in the order encountered, so detailed indices are not required.
		* Before proceeding, the program will reset all molecule 
		indices in the m.r.a structures.  If the atom and residue 
		pointers do not point there, results may be unexpected.
	4. Frees all the molecule space within A.  
	5. Resets the a** and r** pointers to the appropriate locations in E.  
		Locations are based on the molecule indices saved earlier.
	6. Sets a new assembly double-pointer inside E to point to A[0].
	
	This is also why A should not already be a member of E.  If you
	want to duplicate an assembly, this is not the function to use.  
*/

void add_assembly_to_ensemble(
	assembly *A, ///< pointer to the assembly being added (see docs)
	ensemble *E ///< pointer to the ensemble being grown
	){
int i; ///< generic counter
int j; ///< generic counter
int mi; ///< molecule counter
int ri; ///< residue counter
int ai; ///< atom counter
int Enmi; ///< the initial value of nm in E
molindex *rmoli; ///< list of direct indices for A's **r pointers
molindex *amoli; ///< list of direct indices for A's **a pointers


// check sanity of the Assembly for our purposes
// while we're in here, set all molecule indices to what we need them to be 
// If A's information is not arranged in M.r.a format, the function will currently complain and exit.  
if(A[0].nm==0) mywhine("add_assembly_to_ensemble: cannot add assembly with no molecules to an ensemble");
i=j=0;
for(mi=0;mi<A[0].nm;mi++){
	i+=A[0].m[mi][0].nr; // get a count of all residues
	for(ri=0;ri<A[0].m[mi][0].nr;ri++){
// Before proceeding, the program will reset all molecule indices in the m.r.a structures. 
		A[0].m[mi][0].r[ri].moli.i=-1;
		A[0].m[mi][0].r[ri].moli.m=mi;
		A[0].m[mi][0].r[ri].moli.r=ri;
		A[0].m[mi][0].r[ri].moli.a=-1; 
		j+=A[0].m[mi][0].r[ri].na; // get a count of all atoms
		for(ai=0;ai<A[0].m[mi][0].r[ri].na;ai++){ // set molecule indices for the atoms
			A[0].m[mi][0].r[ri].a[ai].moli.i=-1;
			A[0].m[mi][0].r[ri].a[ai].moli.m=mi;
			A[0].m[mi][0].r[ri].a[ai].moli.r=ri;
			A[0].m[mi][0].r[ri].a[ai].moli.a=ai; 
			}
		}
	}
// If A's information is not arranged in M.r.a format, the function will currently complain and exit.  
if(i==0) mywhine("add_assembly_to_ensemble: all molecules in assemnly have no residues; cannot add to ensemble");
if(j==0) mywhine("add_assembly_to_ensemble: all residues in assemnly have no atoms; cannot add to ensemble");

// Make a list of the molecule indices for each atom in the double-pointed atom and residue lists of A. 
if(A[0].nr!=0){
	if(A[0].nr!=i){printf("add_assembly_to_ensemble: raw residue count does not equal A[0].nr (for **r).  Ignoring.\n");}
	if(A[0].nr>i){printf("add_assembly_to_ensemble WARNING: raw residue count is LESS than A[0].nr (for **r).  Ignoring.\n");}
	rmoli=(molindex *)calloc(A[0].nr,sizeof(molindex));
	for(ri=0;ri<A[0].nr;ri++){rmoli[ri]=A[0].r[ri][0].moli;}
	}
if(A[0].na!=0){
	if(A[0].na!=j){printf("add_assembly_to_ensemble: raw atom count does not equal A[0].na (for **a).  Ignoring.\n");}
	if(A[0].na>j){printf("add_assembly_to_ensemble WARNING: raw atom count is LESS than A[0].na (for **a).  Ignoring.\n");}
	amoli=(molindex*)calloc(A[0].na,sizeof(molindex));
	for(ai=0;ai<A[0].na;ai++){amoli[ai]=A[0].a[ai][0].moli;}
	}

// record the offset between A's m.r.a sets and E's. 
Enmi=E[0].nm;
// Allocates and/or reallocates E so that it can accomodate the information found in A.
if(E[0].nm==0) E[0].m=(molecule*)calloc(1,sizeof(molecule)); // initialize if empty
E[0].nm+=A[0].nm; // increment E[0].nm to accomodate new molecules
E[0].m=(molecule*)realloc(E[0].m,E[0].nm*sizeof(molecule)); 

//printf("Enmi is %d; E[0].nm is %d\n",Enmi,E[0].nm);

// Copy atom information in A into the relevant locations in E.
// The molecules are all added in the order encountered, so detailed indices are not required.
for(mi=Enmi;mi<E[0].nm;mi++){
	E[0].m[mi]=A[0].m[mi-Enmi][0]; // copy whole molecule over
	// Free this molecule space within A.  
	free(A[0].m[mi-Enmi]); 
	A[0].m[mi-Enmi]=&E[0].m[mi];
	// reset molecule indicies for residues and atoms
	for(ri=0;ri<E[0].m[mi].nr;ri++){
		E[0].m[mi].r[ri].moli.m+=Enmi;
		for(ai=0;ai<E[0].m[mi].r[ri].na;ai++){ 
			E[0].m[mi].r[ri].a[ai].moli.m+=Enmi;
			}
		}
	}

// Reset the a** and r** pointers to the appropriate locations in E. 
//	It bases the locations on the molecule indices it saved earlier.
for(ri=0;ri<A[0].nr;ri++){ A[0].r[ri]=&E[0].m[rmoli[ri].m].r[rmoli[ri].r]; }
for(ai=0;ai<A[0].na;ai++){ A[0].a[ai]=&E[0].m[amoli[ai].m].r[amoli[ai].r].a[amoli[ai].a]; }

// Set a new assembly double-pointer inside E to point to A[0].
E[0].nA++;
if(E[0].nA==1) E[0].A=(assembly**)calloc(1,sizeof(assembly*)); // initialize if empty
else E[0].A=(assembly**)realloc(E[0].A,E[0].nA*sizeof(assembly*));
E[0].A[E[0].nA-1]=&A[0];

free(rmoli);
free(amoli);

return;
}
