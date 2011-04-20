// Function written by B. Lachele Foley, 2007
#include <mylib.h>
//#include "../inc/mylib.h"

/***************** sscan_file *******************/
/* scans file for number of matches to string s */ 
int sscan_file(FILE *FS,const char *s){
int sfr=0,sfck=1;
char sfs[250];

sfck=fscanf(FS,"%s",sfs);
while (sfck!=EOF) {
	if(strcmp(sfs,s)==0){
		sfr++;
		}
	sfck=fscanf(FS,"%s",sfs);
	}
rewind(FS);

return sfr;
}

/***************** cscan_file *******************/
/* scans file for number of matches to character c*/ 
int cscan_file(FILE *FC,char c){
int cfr=0;
char cfck=0;

cfck=fgetc(FC);
while (cfck!=EOF) {
	if(cfck==c){cfr++;} 
	cfck=fgetc(FC);
	}
rewind(FC);

return cfr;
}
