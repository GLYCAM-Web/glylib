/************************** init_struct() ******************************/
/* This is load_pdb's main initialization function.  It assigns default rule
 * values to some data structures and calls other functions to assign
 * values to other functions.
 * Author: Lachele Foley
 * Modified for use in library by: Michael Tessier
 */
#include <load_pdb.h>
//#include "../inc/load_pdb.h"
extern int DEBUG,INWC;
void init_struct(){
int isa=0, isb=0,isc=0, iswc=0;
char istst='\0';
if(DEBUG>=0){printf("init_struct at 1 \n");} 
// Get number of lines in input file
rewind(IN);
while(istst!=EOF){
	istst=fgetc(IN);
	if(istst=='\n') iswc++;
	}
INWC=iswc;

/* Allocate array for line info */
ln=(linedef*)calloc(iswc,sizeof(linedef)); 
for(isa=0;isa<iswc;isa++){
	ln[isa].m='\0'; 
	ln[isa].n=isa+1;
	ln[isa].ignore=1;
	for(isb=0;isb<25;isb++){
		for(isc=0;isc<81;isc++){
			ln[isa].f[isb].c[isc]='\0';
			}
		}
	}
if(DEBUG>=0){printf("init_struct at  9\n");}
/* Read in values defining line and field types */
pdb_def();

if(DEBUG>=0){printf("init_struct at 12 \n");}
return;
}


