/** \file get_moiety_selection_assembly.c Returns an assembly pointing to the selections 
	indicated by the entries in a set of moiety_selections.
	
	If two different moiety_selections overlap (select some of the same items),
	the assembly will not contain redundant entries.  For example, if residue A
	is selected twice in one moiety selection or if it is selected in more than
	one moiety, the assembly returned will contain one pointer to residue A.  If
	choosing by name, and multiple residues are named "A", there will only be one
	pointer to each residue of that name no matter how many times "A" is specified
	as a selection.

	Note that a negative moiety_selection does not unselect -- it merely specifies
	a selection by saying "select all except".  Use the function (not written yet)
	called remove_moiety_selection_from_assembly to unselect items.

	The assembly will consist only of pointers to any members any classes represented 
	in the moiety selection.  For example, if a moiety selection contains residues
	A, B and C, there will be pointers to A, B and C, but not to any of the atoms
	contained within those residues.  However, if the moiety selections also contain
	explicit mention of atoms from any of those residues, the assembly will also
	contain pointers to those atoms explicitly mentioned.  It is up to some other
	function to prune or clean up or expand, etc.

	NOTES:
	- This function will set all the moli's for the residues and atoms and the mi
		member in the molecule structure.  It will not check to see if you had
		something else there first.  If you need to keep some information other 
		than the standard indices, find somewhere else to put it.
	- The entries in the assembly will be ordered as the order in the original
		ensemble.  For example, if the moiety_selection chooses residues number
		5, 250, 48, 2 and 143, the assembly will contain pointers to residues
		number 2, 5, 48, 143 and 250, in that order.  
	- If residues spanning different molecules are selected, the assembly will
		not group them by molecule.  They will be in one long list of 
		residue pointers.  You can use the moli entries to determine membership
		in a particular molecule, if needed.
	- Of course, name, number or index entries in MS must exactly match names,
		numbers and/or indices of entries in E.  If not, however, the function
		will not complain.  It simply will not add selections to the assembly
		where it does not find a match.
	- Of course, the moiety_selection should not contain any uninitialized or
		members.  
*/

#include "mylib.h"
#include "molecules.h"


//START HERE -- change name to get_moiety_selection_assembly_from_ensemble
// Make similar called get_moiety_selection_assembly_from_assembly
assembly get_moiety_selection_assembly(ensemble *E, int nMS, moiety_selection *MS){
int 	am=0,ar=0,aa=0,ai=0,ams=0, ///< counters
	aflag=0, ///< Flag used for setting negative selections
	numM=0, ///< Number of molecule pointers needed for assembly
	numR=0, ///< Number of residue pointers needed for assembly
	numA=0; ///< Number of atom pointers needed for assembly
assembly A; ///< the assembly to return
ensemble_tree_index IA; ///< the indices for holding selections

// Set up the tree of indices (IA) analogous to the m,r,a structure in E
//	initialize them all to be 0
IA.nm=E[0].nm;
IA.mi=(int*)calloc(IA.nm,sizeof(int));
IA.m=(molecule_tree_index*)calloc(IA.nm,sizeof(molecule_tree_index));
for(am=0;am<IA.nm;am++){ // for each molecule
	IA.m[am].nr=E[0].m[am].nr;
	IA.m[am].ri=(int*)calloc(IA.m[am].nr,sizeof(int));
	IA.m[am].r=(residue_tree_index*)calloc(IA.m[am].nr,sizeof(residue_tree_index));
	for(ar=0;ar<IA.m[am].nr;ar++){ // for each residue 
		IA.m[am].r[ar].na=E[0].m[am].r[ar].na;
		IA.m[am].r[ar].ai=(int*)calloc(IA.m[am].r[ar].na,sizeof(int));
		}
	}


/** Walk through each m, r, a of the ensemble.  Look for matches.  
*
* If there is a match, change the corresponding index in IA to equal one.  
* Don't check to see whether it is already one -- that doesn't matter. */

for(ams=0;ams<nMS;ams++){

/// For each positive moiety_selection passed to the function
if(MS[ams].posneg==1){ 
	// Molecules
	// For each molecule index specified
	for(ai=0;ai<MS[ams].nmi;ai++){ IA.mi[MS[ams].mi[ai]]=1; }
	// Check each molecule for names and numbers
	for(am=0;am<IA.nm;am++){ 
		// For each name specified
		for(ai=0;ai<MS[ams].nmN;ai++){ if((E[0].m[am].N!=NULL)&&(strcmp(E[0].m[am].N,MS[ams].mN[ai])==0)){IA.mi[am]=1;} }
		// For each number specified
		for(ai=0;ai<MS[ams].nmn;ai++){ if(E[0].m[am].n==MS[ams].mn[ai]){IA.mi[am]=1;} }	
	
		// Residues
		// For each index specified
		for(ai=0;ai<MS[ams].nri;ai++){ IA.m[am].ri[MS[ams].ri[ai]]=1; }
		// Now, check each residue for names and numbers 
		for(ar=0;ar<IA.m[am].nr;ar++){ 
			// For each name specified
			for(ai=0;ai<MS[ams].nrN;ai++){ if((E[0].m[am].r[ar].N!=NULL)&&(strcmp(E[0].m[am].r[ar].N,MS[ams].rN[ai])==0)){IA.m[am].ri[ar]=1;} }
			// For each number specified
			for(ai=0;ai<MS[ams].nrn;ai++){ if(E[0].m[am].r[ar].n==MS[ams].rn[ai]){IA.m[am].ri[ar]=1;} }	

			// Atoms
			// For each index specified
			for(ai=0;ai<MS[ams].nai;ai++){ IA.m[am].r[ar].ai[MS[ams].ai[ai]]=1; }
			// Now, check each residue for names and numbers 
			for(aa=0;aa<IA.m[am].r[ar].na;aa++){ 
// For each atom name specified
for(ai=0;ai<MS[ams].naN;ai++){ if((E[0].m[am].r[ar].a[aa].N!=NULL)&&(strcmp(E[0].m[am].r[ar].a[aa].N,MS[ams].aN[ai])==0)){IA.m[am].r[ar].ai[aa]=1;} }
// For each atom number specified
for(ai=0;ai<MS[ams].nan;ai++){ if(E[0].m[am].r[ar].a[aa].n==MS[ams].an[ai]){IA.m[am].r[ar].ai[aa]=1;} } 
				} // close aa loop
			} // close ar loop
		} // close am loop
	} // close if positive MS condition

/// For each negative moiety_selection passed to the function
if(MS[ams].posneg==1){ 
	// Molecules
	// Check each molecule for names and numbers
	for(am=0;am<IA.nm;am++){ 
		// For each molecule index specified
		aflag=1;
		for(ai=0;ai<MS[ams].nmi;ai++){ if(am==MS[ams].mi[ai]) aflag=0; }
		if(aflag==1) { IA.mi[am]=1; }
		// For each name specified
		aflag=1;
		for(ai=0;ai<MS[ams].nmN;ai++){ if((E[0].m[am].N!=NULL)&&(strcmp(E[0].m[am].N,MS[ams].mN[ai])==0)){aflag=0;} }
		if(aflag==1) { IA.mi[am]=1; }
		// For each number specified
		aflag=1;
		for(ai=0;ai<MS[ams].nmn;ai++){ if(E[0].m[am].n==MS[ams].mn[ai]){aflag=0;} }	
		if(aflag==1) { IA.mi[am]=1; }
	
		// Residues
		// Now, check each residue for names and numbers 
		for(ar=0;ar<IA.m[am].nr;ar++){ 
			// For each index specified
			aflag=1;
			for(ai=0;ai<MS[ams].nri;ai++){ if(ar==MS[ams].ri[ai]) aflag=0; }
			if(aflag==1) { IA.m[am].ri[ar]=1; }
			// For each name specified
			aflag=1;
			for(ai=0;ai<MS[ams].nrN;ai++){ if((E[0].m[am].r[ar].N!=NULL)&&(strcmp(E[0].m[am].r[ar].N,MS[ams].rN[ai])==0)){aflag=0;} }
			if(aflag==1) { IA.m[am].ri[ar]=1; }
			// For each number specified
			aflag=1;
			for(ai=0;ai<MS[ams].nrn;ai++){ if(E[0].m[am].r[ar].n==MS[ams].rn[ai]){aflag=0;} }	
			if(aflag==1) { IA.m[am].ri[ar]=1; }

			// Atoms
			// Now, check each residue for names and numbers 
			for(aa=0;aa<IA.m[am].r[ar].na;aa++){ 
// For each index specified
aflag=1;
for(ai=0;ai<MS[ams].nai;ai++){ if(aa==MS[ams].ai[ai]) aflag=0; }
if(aflag==1) { IA.m[am].r[ar].ai[aa]=1; } 
// For each atom name specified
aflag=1;
for(ai=0;ai<MS[ams].naN;ai++){ if((E[0].m[am].r[ar].a[aa].N!=NULL)&&(strcmp(E[0].m[am].r[ar].a[aa].N,MS[ams].aN[ai])==0)){aflag=0;} }
if(aflag==1) { IA.m[am].r[ar].ai[aa]=1; } 
// For each atom number specified
aflag=1;
for(ai=0;ai<MS[ams].nan;ai++){ if(E[0].m[am].r[ar].a[aa].n==MS[ams].an[ai]){aflag=0;} } 
if(aflag==1) { IA.m[am].r[ar].ai[aa]=1; } 
				} // close aa loop
			} // close ar loop
		} // close am loop
	} // close if negative MS condition

	} // close ams loop

// after all the indices are set for selection, count, in IA, the number of
// 	m, r and a pointers the assembly will need and allocate them.
numM=numR=numA=0;
for(am=0;am<IA.nm;am++){ if(IA.mi[am]>0) numM++;
	for(ar=0;ar<IA.m[am].nr;ar++){ if(IA.m[am].ri[ar]>0) numR++;
		for(aa=0;aa<IA.m[am].r[ar].na;am++){ if(IA.m[am].r[ar].ai[aa]==1) numA++;
			}
		}
	}
A.nm=numM;
A.m=(molecule**)calloc(A.nm,sizeof(molecule*));
A.nr=numR;
A.r=(residue**)calloc(A.nr,sizeof(residue*)); 
A.na=numA;
A.a=(atom**)calloc(A.na,sizeof(atom*)); 

// Walk through IA and set an appropriate pointer whenever a "1" is encountered
numM=numR=numA=0;
for(am=0;am<IA.nm;am++){ 
	if(IA.mi[am]>0){
		A.m[numM]=&E[0].m[am];
		numM++;
		if(numM>A.nm){mywhine("get_moiety_selection_assembly: memory allocation issue in A.nm/numM");}
		}
	for(ar=0;ar<IA.m[am].nr;ar++){ 
		if(IA.m[am].ri[ar]>0){
			A.r[numR]=&E[0].m[am].r[ar];
			numR++;
			if(numR>A.nr){mywhine("get_moiety_selection_assembly: memory allocation issue in A.nr/numR");}
			}
		for(aa=0;aa<IA.m[am].r[ar].na;am++){ 
			if(IA.m[am].r[ar].ai[aa]==1){
				A.a[numA]=&E[0].m[am].r[ar].a[aa];
				numA++;
				if(numA>A.na){mywhine("get_moiety_selection_assembly: memory allocation issue in A.na/numA");}
				}
			}
		}
	}

return A;
}
