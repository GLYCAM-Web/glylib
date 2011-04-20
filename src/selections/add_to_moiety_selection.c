/** \file add_to_moiety_selection.c  Parses input string and adds selection to 
	a "moiety_selection" structure.

This function might not be complete -- it is being written for a certain project,
and should be expanded as needed.

NOTE!!!  The moiety_selection *MUST* be allocated and initialized before calling 
	this function.  Initialization means that sub-arrays must have had an 
	initial allocation (so that realloc won't choke).

	It also must be initialized so that the "nX" variables start at zero. 
	They need not be zero when they come to this function, but they need 
	to be something real and not whatever garbage was left in memory.


The current syntax is:

	number  moiety-size-class  moiety-designation  optional_NOT  designation_list

Only one such 'sentence' is allowed per call.

Syntax Examples:

Add residues number 1, 2, 8 and 38 to the selection:

	4 Residue Numbers 1 2 8 38

Add atoms not named O5 and O6 

	2 Atoms Names -NOT- O5 O6  

NOTE:  Regarding entries, the function inspects only the smallest number of characters
	necessary to differentiate keywords.  So, for the moiety-size-class entry, the
	strings "Residue" "R" "Raisin" and "RESI" are all the same to the program.

Acceptable values for input.  The star (*) indicates the part of the function that will 
	be written now (starting 20080910-ish):

	Number:  integer
		Eventually, this could be written to also read "range", but not today.

	moiety-size-class:	Atom  Residue(*)  Molecule  Assembly Ensemble Group
		First letter differentiates.

	moiety-designation: 	Name  Number  Index
		First two letters differentiate (actually, using second 
			character in the string presently).
		Number = the number in the input data set (e.g., pdb file).
		Index = the array index as stored internally to the program

	optional_NOT:	-NOT- 
		Must include all five characters, case sensitive.
		A single moiety_selection must be all positive or all negative.
		The calling function should decide before adding info to it.
		If SEL doesn't match, an error will be produced.

	designation_list:  (varies)
		Space-separated list of names or numbers.  Names must exactly 
		match the string contained in "N".  Currently, spaces are not 
		permitted in the names. (someone could re-write to allow it...)

More about the selections:  
	- Selections will be set only within the relevant class.  For example, 
		if the class is set to "residue", then only residues will be 
		added to the structure; their associated atoms will not be 
		explicitly included.  Similarly, if the -NOT- option is chosen 
		with class residue, all residues except those listed will be 
		added to the selection, but not the associated atoms.
	- Only the first "number" designations will be scanned from the
		designation_list.  Any trailing characters or words will be
		ignored after the first "number" words are found.  
	- Empty selection strings will be ingored and no error returned.  The
		moiety_selection structure will be unaltered.
	- Errors will be returned if:
		- The first word in SEL is not an integer
		- Fewer than "number" words are found in the designation list
		- No match is found for either class or designation
		- The existing moiety_selection is set to be a positive list
			but SEL is negative or vice-versa.
		- Note, of course, that this function merely parses a 
			selection string.  It has no knowledge whether the
			molecule, etc., actually contains the selected items.

*/

#include <mylib.h>
#include <molecules.h>

void add_to_moiety_selection(moiety_selection *M, const char *SEL){
int 	aa=0, ///< counter
	ab=0, ///< counter
	ac=0, ///< counter
	posneg=0, ///< flag for negative selection
	number=0; ///< number of designations in the list
char 	*tp, ///< Temporary pointer
	*ts, ///< Temporary string
	*S, ///< The trimmed selection
	*class, ///< The class
	*des, ///< The designation
	**d; ///< The designation_list

/// Make sure there is no leading or trailing whitespace in S to cause confusion
S=strdup(prune_string_whitespace(SEL));

/// grab three items, place into number, class and designation
tp=S; // use a temporary pointer (tp) to move around in S.

// Get Integer
aa=strcspn(tp," \t\n\v\r\f"); ///< Find length of first word -- this should be the integer
if(aa==0){return;} ///< If the string S is empty, return without complaint
ts=(char*)calloc(aa+1,sizeof(char));
strncpy(ts,tp,aa); 
ts[aa]='\0';
ab=sscanf(ts,"%d",&number); ///< Scan in the first integer
if(ab==0){mywhine("add_to_moiety_selection: Initial item in SEL not an integer");}
free(ts);

// Get Class
tp+=aa; // move the pointer past the first item
aa=strspn(tp," \t\n\v\r\f");
tp+=aa; // move past the whitespace
aa=strcspn(tp," \t\n\v\r\f"); ///< Find length of second word -- this should be the class
class=(char*)calloc((aa+1),sizeof(char));
strncpy(class,tp,aa); 
class[aa]='\0';



// Get Designation
tp+=aa; // move the pointer past the second item
aa=strspn(tp," \t\n\v\r\f");
tp+=aa; // move past the whitespace
aa=strcspn(tp," \t\n\v\r\f"); ///< Find length of third word -- this should be the designation
des=(char*)calloc((aa+1),sizeof(char));
strncpy(des,tp,aa); 
des[aa]='\0';

// See if there is an Optional -NOT-
ts=strstr(tp,"-NOT-"); 
if(ts!=NULL){ // if a not entry was found
	posneg=-1;
	tp=ts+5; // move the tp pointer past the "-NOT-"
	}
else{posneg=+1;}

if((M[0].posneg!=0)&&(M[0].posneg!=posneg)){mywhine("add_to_moiety_selection: positive/negative mismatch with SEL and M.");}
M[0].posneg=posneg; // In case it was zero

// separate the designation_list by scanning for non-whitespace
d=(char**)calloc(number,sizeof(char*));
for(ab=0;ab<number;ab++){
	aa=strspn(tp," \t\n\v\r\f");
	tp+=aa; // move past the whitespace
	aa=strcspn(tp," \t\n\v\r\f"); ///< Find length of the ab'th item in the designation_list
	d[ab]=(char*)calloc((aa+1),sizeof(char));
	strncpy(d[ab],tp,aa); 
	d[ab][aa]='\0'; 
	}

// add the string entries to the appropriate parts of the moiety_selection
switch(class[0]){ ///< First letter of class differentiates molecule, residue or atom
	case 'M':
	case 'm':
		switch(des[1]){ ///< Second letter of designation differentiates name, number or index
			case 'A':
			case 'a':
				ab=M[0].nmN;
				M[0].nmN+=number;
				M[0].mN=(char**)realloc(M[0].mN,M[0].nmN*sizeof(char*));
				for(aa=0;aa<number;aa++){M[0].mN[aa+ab]=strdup(d[aa]);}
				break;
			case 'U':
			case 'u':
				ab=M[0].nmn;
				M[0].nmn+=number;
				M[0].mn=(int*)realloc(M[0].mn,M[0].nmn*sizeof(int));
				for(aa=0;aa<number;aa++){ac=sscanf(d[aa],"%d",&M[0].mn[aa+ab]);
					if(ac!=1){mywhine("add_to_moiety_selection: expected integer for des-number, got something else");}}
				break;
			case 'N':
			case 'n':
				ab=M[0].nmi;
				M[0].nmi+=number;
				M[0].mi=(int*)realloc(M[0].mi,M[0].nmi*sizeof(int));
				for(aa=0;aa<number;aa++){ac=sscanf(d[aa],"%d",&M[0].mi[aa+ab]);
					if(ac!=1){mywhine("add_to_moiety_selection: expected integer for index, got something else");}}
				break;
			default:
				mywhine("add_to_moiety_selection: unrecognized moiety-designation");
			}
		break;
	case 'R':
	case 'r': 
		switch(des[1]){ ///< Second letter of designation differentiates name, number or index
			case 'A':
			case 'a':
				ab=M[0].nrN;
				M[0].nrN+=number;
				M[0].rN=(char**)realloc(M[0].rN,M[0].nrN*sizeof(char*));
				for(aa=0;aa<number;aa++){M[0].rN[aa+ab]=strdup(d[aa]);}
				break;
			case 'U':
			case 'u':
				ab=M[0].nrn;
				M[0].nrn+=number;
				M[0].rn=(int*)realloc(M[0].rn,M[0].nrn*sizeof(int));
				for(aa=0;aa<number;aa++){ac=sscanf(d[aa],"%d",&M[0].rn[aa+ab]);
					if(ac!=1){mywhine("add_to_moiety_selection: expected integer for des-number, got something else");}}
				break;
			case 'N':
			case 'n':
				ab=M[0].nri;
				M[0].nri+=number;
				M[0].ri=(int*)realloc(M[0].ri,M[0].nri*sizeof(int));
				for(aa=0;aa<number;aa++){ac=sscanf(d[aa],"%d",&M[0].ri[aa+ab]);
					if(ac!=1){mywhine("add_to_moiety_selection: expected integer for index, got something else");}}
				break;
			default:
				mywhine("add_to_moiety_selection: unrecognized moiety-designation");
			}
		break;
	case 'A':
	case 'a': 
		switch(des[1]){ ///< Second letter of designation differentiates name, number or index
			case 'A':
			case 'a':
				ab=M[0].naN;
				M[0].naN+=number;
				M[0].aN=(char**)realloc(M[0].aN,M[0].naN*sizeof(char*));
				for(aa=0;aa<number;aa++){M[0].aN[aa+ab]=strdup(d[aa]);}
				break;
			case 'U':
			case 'u':
				ab=M[0].nan;
				M[0].nan+=number;
				M[0].an=(int*)realloc(M[0].an,M[0].nan*sizeof(int));
				for(aa=0;aa<number;aa++){ac=sscanf(d[aa],"%d",&M[0].an[aa+ab]);
					if(ac!=1){mywhine("add_to_moiety_selection: expected integer for des-number, got something else");}}
				break;
			case 'N':
			case 'n':
				ab=M[0].nai;
				M[0].nai+=number;
				M[0].ai=(int*)realloc(M[0].ai,M[0].nai*sizeof(int));
				for(aa=0;aa<number;aa++){ac=sscanf(d[aa],"%d",&M[0].ai[aa+ab]);
					if(ac!=1){mywhine("add_to_moiety_selection: expected integer for index, got something else");}}
				break;
			default:
				mywhine("add_to_moiety_selection: unrecognized moiety-designation");
			}
		break;
	default:
		mywhine("add_to_moiety_selection: unrecognized moiety-size-class");
	}

free(ts);
free(class);
free(des);
for(aa=0;aa<number;aa++){free(d[aa]);}
free(d);

return;
}
