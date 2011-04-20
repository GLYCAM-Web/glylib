/* \file fileslurp.c 
\addtogroup FILE_UTILS
\brief	Read an entire file into a structure to be processed
	later within the calling function.

	Argument: A fileset structure.

	Returns: Fileslurp structure containing each line in
		the file (with newline intact).

	Begun on 20080706 by BLFoley.
*/

#include <mylib.h>
#include <general.h>
//fileslurp slurp_file(fileset F);


fileslurp slurp_file(fileset F){
fileslurp S;
char templine[2001],*line;
int doneflag=0;

// fgets each line in the file, place in temporary hold
// issue warning if any line is longer than 2000 chars, but
// finish reading the thing in anyhow.

// open the file
F.F=myfopen(F.N,"r");
// read the file
S.n=0;
S.L=(char**)calloc(1,sizeof(char*));
while(fgets(templine,2000,F.F)!=NULL){
	line=strdup(templine);
	if(strchr(templine,'\n')==NULL){ // got a line over 2000 chars
		fprintf(stderr,"found line over 2000 chars in file ");
		fprintf(stderr,"%s at line %d. reading anyway.\n",F.N,(S.n+1));
		doneflag=1;
		while(doneflag==1){
			fgets(templine,2000,F.F);
	line=(char*)realloc(line,(strlen(line)+strlen(templine)+1)*sizeof(char));
			strcat(line,templine);
			if(strchr(templine,'\n')!=NULL){doneflag=0;}
			}
		}
	S.n++;
	S.L=(char**)realloc(S.L,S.n*sizeof(char*));
	S.L[S.n-1]=strdup(line);
	} // close the while loop

return S;
} 
