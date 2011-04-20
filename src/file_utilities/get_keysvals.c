/** \file get_keysvals.c Functions for extracting lists of keywords 
and values.

Contains functions to ease reading of files (slurped or not).
Begun on 20080816 by BLFoley.
*/

#include "mylib.h"
//#include "general.h"

/****************** get_keysvals_from_slurp() ******************/
/** Reads in keyword/value pairs from a file, drops them in a
structure designed for the purpose.

Some notes:

	- If ignore_hash is 0, everything from a hash mark (#) to the end of the line is ignored.
	- Any line not containing the separator SEP is ignored.
	- Keywords and values may not be separated by a newline.
	- If newline_sep is set to zero, everything* from the start of the line to the first instance 
		of SEP will be regarded as the KEYWORD and everything* from the first instance of SEP
		to the end of the line will be regarded as the VALUE. *Intial and final whitespace will 
		be removed from both keyword and value, but all else, including internal whitespace,
		will be retained.
	- If newline_sep is not set to zero, neither keywords nor values may contain whitespace.
	- SEP should be passed as a string complete with a null terminator.
*/
gly_keysvals get_keysvals_from_slurp(fileslurp F, const char *SEP, int newline_sep, int ignore_hash){
int ga=0; ///< counter
char *t; ///< temporary pointer
gly_keysvals KV;

KV.n=0;
KV.K=(char**)calloc(1,sizeof(char*));
KV.V=(char**)calloc(1,sizeof(char*));

for(ga=0;ga<F.n;ga++){ ///< Loop through the slurp
	if(ignore_hash==0){
		t=strchr(F.L[ga],'#'); ///< See if there is a hash in this line
		if(t!=NULL) t[0]='\0'; ///< If there is a hash, call that the end of the string
		}
	if(strstr(F.L[ga],SEP)==NULL) continue; ///< If no separator this line, ignore it
	if(newline_sep==0){///< If newline_sep is set
		KV.n++;
		KV.K=(char**)realloc(KV.K,KV.n*sizeof(char*));
		KV.V=(char**)realloc(KV.V,KV.n*sizeof(char*));
		t=strstr(F.L[ga],SEP); ///< Split string at SEP (find with strstr)
		t[0]='\0'; ///< Set end of the keyword string
		KV.K[KV.n-1]=strdup(prune_string_whitespace(F.L[ga])); ///< Prune whitespace and copy to K
		t+=(strlen(SEP)); ///< Move pointer to the end of the separator
		KV.V[KV.n-1]=strdup(prune_string_whitespace(t)); ///< Prune whitespace and copy to V
		}
	else{
		fprintf(stderr,"\n\nIn get_keysvals_from_slurp (or from file):\n");
		fprintf(stderr,"The use of multiple keyword/value pairs on one line is not yet supported.\n");
		fprintf(stderr,"Feel free to write it in...   Exiting now.\n\n");
		exit(0);
/** 
Else
	- Grab K
	- Check for SEP 
	- Grab V
	--- Allow for line-crossing?  No!!! (what a pain...)
*/
		}
	} // close loop through the slurp


return KV;
}
/****************** get_keysvals_from_file() ******************/
/** Reads in keyword/value pairs from a file, drops them in a
structure designed for the purpose.
*/
gly_keysvals get_keysvals_from_file(fileset F, const char *SEP, int newline_sep, int ignore_hash){
fileslurp FS;
gly_keysvals KV;

FS=slurp_file(F); ///< Slurp file 
KV=get_keysvals_from_slurp(FS, SEP, newline_sep, ignore_hash); ///< Call the slurp version...

return KV;
}
