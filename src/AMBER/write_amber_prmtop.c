/** \file write_amber_prmtop.c Writes out an amber prmtop file.

There should be an amber_prmtop structure available that contains the proper
information.  It is up to the calling function to make any needed modifications
to that structure.  See also the read_amber_prmtop function.  The typical flow
goes something like this:

 -# Read an amber prmtop file using read_amber_prmtop
 	- The entire contents of the original file can be found in the void pointer.
		See the documentation on use of the void pointer.
 -# Decide which parts to change and which parts of the original file to include
	in the output file.
	- Note that not all parts of the original file are necessarily parsed into
		the assembly structure.  See parse_amber_prmtop for more info.
 -# Copy the updated prmtop contents to an amber_prmtop struct.  
	- It is probably safest to make a new amber-prmtop struct for the changed info.
	- If memory space is an issue, the information can be written back into the 
		void pointer.  Alternately, and better, a new prmtop can be generated 
		and the old one (in the void pointer) can be freed.
 -# Call this function to write an amber prmtop file out.

Note, also, that the amber_prmtop struct contains content and formatting information.
It is the duty of the calling program to get that right.  This function will provide
only internal checks -- that it, it will make sure that the string lengths match the
data information.  **But** it will get formatting information from the interpreted
units in the amber_prmtop_section structure and -not- from the original formatting
string.  **However** it will write the original string in *FORMAT out to the file. So,
make sure there is a match.
*/

#include <AMBER/amber.h>

void write_amber_prmtop(fileset F, amber_prmtop *P){
int wa=0,wb=0,wc=0,NCL=0;

/// Open the file
F.F=myfopen(F.N,"w");

/// Write the version info
fprintf(F.F,"%%VERSION %s\n",P[0].VERSION);

/** Loop through all the sections 
	
	- Check that each string length matches the format length
	- Check that the number of units matches the length of the array **D
	- Write the required number per lie
	- Write the contents
*/
for(wa=0;wa<P[0].nS;wa++){ // for each section found in the file
	fprintf(F.F,"%%FLAG %s\n",P[0].S[wa].N);
	fprintf(F.F,"%%FORMAT(%s)\n",P[0].S[wa].FORMAT);
	NCL=floor(P[0].S[wa].nt/P[0].S[wa].npl); // the number of full (complete) lines
	wc=0;
	for(wb=0;wb<NCL;wb++){ // write the complete lines
		for(wc=wc;wc<P[0].S[wa].npl;wc++){ // each entry
			if(strlen(P[0].S[wa].D[wc])!=P[0].S[wa].nc){mywhine("Improper format found in write_amber_prmtop");}
			fprintf(F.F,"%s",P[0].S[wa].D[wc]);
			}
		fprintf(F.F,"\n");
		}
	for(wc=wc;wc<P[0].S[wa].nt;wc++){ // write the last incomplete line
		if(strlen(P[0].S[wa].D[wc])!=P[0].S[wa].nc){mywhine("Improper format found in write_amber_prmtop");}
		fprintf(F.F,"%s",P[0].S[wa].D[wc]); 
		}
	fprintf(F.F,"\n");
	} // close the loop over each section 
fclose(F.F);
return;	
}
