/**************************  rwm_line() *******************************/
/* Author: Lachele Foley
   Modified for use in library by: Michael Tessier*/
#include <load_pdb.h>
//#include "../inc/load_pdb.h"
extern int DEBUG,WARN,SILENT;
/* This reads in the line */
void rwm_line(int rwmln){
int rwmdone=0, rwma=0, rwmb=0, rwmstop=-1, rwmmt=-1, rwmnl=-1;
int shortline=-1;
linetype rwmlt;
//char rwmaCL='I',rwmaRN='I';
char rwmtest='\0';
/* Read in first six characters and determine line type */
if(DEBUG>=1){printf("rwm_line 1\n");}
for(rwma=0;rwma<6;rwma++){
	if(rwmstop!=0){
		rwmdone=fgetc(IN);
if(DEBUG>=2){printf("Just got new rwmdone \n");}
		if(rwmdone==EOF){read_neek("EOF",rwmln,rwma);}
		}
	if(((rwmdone=='\n')||(rwmdone==' '))&&(rwma==0)){
if(DEBUG>=2){printf("Empty line found \n");}
		if(rwmdone=='\n'){rwmnl=0;}
		rwmmt=0;
		}
	if((rwmdone=='\n')||(rwmstop==0)){
if(DEBUG>=2){printf("resetting rwmdone \n");}
		rwmnl=0;
		rwmdone=' ';
		rwmstop=0;
		}
	ln[rwmln-1].f[0].c[rwma]=rwmdone;
	}
if(rwmmt==0){ /* If this is an empty line */
	if(rwmnl!=0){
		rwmdone=fgetc(IN);
		while(rwmdone!='\n'){
			rwmdone=fgetc(IN);
			}
		}
	/* if(WARN!=0){fprintf(OUT,"\n");} */
	return;
	}
if(DEBUG>=1){printf("rwm_line 2\n");}
// if the first six chars seem to be ok...
// null-terminate the string
ln[rwmln-1].f[0].c[6]='\0';

// made it here?  find out the line type
if(DEBUG>=1){printf("rwm_line 3 fieldname is >>>%s<<<\n",ln[rwmln-1].f[0].c);}
rwmlt=get_type(ln[rwmln-1].f[0].c); 
ln[rwmln-1].a=rwmlt.a; // class 
ln[rwmln-1].b=rwmlt.b; // number for the sub-class



/* Assign action switch from internal defaults */
if(DEBUG>=1){printf("rwm_line 4\n");}
if(DEBUG>=1){printf("rwmlt.a is %d and rwmlt.b is %d\n",rwmlt.a,rwmlt.b);} 
/* If get_type returned a valid class: */


//if(rwmlt.a==-1){ // If this line type is unknown
// do something....
	//rwmaCL='I';
	//rwmaRN='D';
	//}
//If there line was less than 6 characters
if(rwmstop == 0)
 return;

// If the rule is to modify, read in the rest of the line, first...  
for(rwma=1;rwma<pdb_a[rwmlt.a].b[rwmlt.b].f;rwma++){
//printf("in the read rest of line loop, chars is %d\n",pdb_a[rwmlt.a].b[rwmlt.b].c[rwma]);
        for(rwmb=0;rwmb<pdb_a[rwmlt.a].b[rwmlt.b].c[rwma];rwmb++){
                ln[rwmln-1].f[rwma].c[rwmb]=fgetc(IN);
//printf("the current character (%d) is >>%c<<\n",rwmb,ln[rwmln-1].f[rwma].c[rwmb]);
                //fprintf(OUTC,"%c",ln[rwmln-1].f[rwma].c[rwmb]);
                if(ln[rwmln-1].f[rwma].c[rwmb]=='\n'){
                        rwmb=pdb_a[rwmlt.a].b[rwmlt.b].c[rwma];
                        rwma=pdb_a[rwmlt.a].b[rwmlt.b].f;
			shortline=0;
			rwmtest='\n';
                        }
                if(ln[rwmln-1].f[rwma].c[rwmb]==EOF){
                        read_neek("EOF",(rwmln),(rwma+1));
                        }
                }
        ln[rwmln-1].f[rwma].c[rwmb]='\0';
        } 

//printf("rwmtest is >>%c<<\t",rwmtest);
if(shortline!=0) rwmtest=fgetc(IN);
//printf("rwmtest is >>%c<<\n",rwmtest);

if((rwmtest!='\n')&&(rwmtest!=EOF)){
	printf("Found line longer than 80 characters (likely line %d).  Ignoring rest of line.\n",rwmln);
	rwmtest='\0';
	while(rwmtest!='\n'){rwmtest=fgetc(IN);}
	}
if(DEBUG>=0){printf("rwm_line 13\n");}
/*if((WARN==0)&&(SILENT!=0)){fflush(OUTC);} */
return;
}

/*This function written by Michael Tessier 20080828*/
void rwm_line_char(char* curLine, int rwmln){
int rwmdone=0, rwma=0, rwmb=0, rwmstop=-1, rwmmt=-1;
int itr,ishorta=0,ishortb=0;
linetype rwmlt;
//char rwmaCL='I',rwmaRN='I';
/* Read in first six characters and determine line type */
if(DEBUG>=1){printf("rwm_line 1\n");}
for(rwma=0;rwma<6;rwma++){
	if(rwmstop!=0){
		rwmdone=curLine[rwma];
if(DEBUG>=2){printf("Just got new rwmdone \n");}
		//if(rwmdone==EOF){read_neek("EOF",rwmln,rwma);}
		}
	if(((rwmdone=='\n')||(rwmdone==' '))&&(rwma==0)){
if(DEBUG>=1){printf("Empty line found \n");}
		rwmmt=0;
		}
	if((rwmdone=='\n')||(rwmstop==0)){
if(DEBUG>=2){printf("resetting rwmdone \n");}
		rwmdone=' ';
		rwmstop=0;
		}
	ln[rwmln-1].f[0].c[rwma]=rwmdone;
	}
if(rwmmt==0){ /* If this is an empty line */
	/*if(WARN!=0){fprintf(OUT,"\n");}*/
	return;
	}
if(DEBUG>=1){printf("rwm_line 2\n");}
// if the first six chars seem to be ok...
// null-terminate the string
ln[rwmln-1].f[0].c[6]='\0';

// made it here?  find out the line type
if(DEBUG>=1){printf("rwm_line 3 fieldname is >>>%s<<<\n",ln[rwmln-1].f[0].c);}
rwmlt=get_type(ln[rwmln-1].f[0].c); 
ln[rwmln-1].a=rwmlt.a; // class 
ln[rwmln-1].b=rwmlt.b; // number for the sub-class



if(DEBUG>=1){printf("rwmlt.a is %d and rwmlt.b is %d\n",rwmlt.a, rwmlt.b);}
/* Assign action switch from internal defaults */
if(DEBUG>=1){printf("rwm_line 4\n");}
//if(DEBUG>=0){printf("rwmlt.a is %d and rwmlt.b is %d\n",rwmlt.a,rwmlt.b);} 
/* If get_type returned a valid class: */


//if(rwmlt.a==-1){ // If this line type is unknown
// do something....
	//rwmaCL='I';
	//rwmaRN='D';
	//}
//If there line was less than 6 characters
if(rwmstop == 0)
 return;
itr = 6;
// If the rule is to modify, read in the rest of the line, first...  
for(rwma=1;rwma<pdb_a[rwmlt.a].b[rwmlt.b].f;rwma++){
        for(rwmb=0;rwmb<pdb_a[rwmlt.a].b[rwmlt.b].c[rwma];rwmb++){
                ln[rwmln-1].f[rwma].c[rwmb]=curLine[itr];
                //fprintf(OUTC,"%c",ln[rwmln-1].f[rwma].c[rwmb]);
                if(ln[rwmln-1].f[rwma].c[rwmb]=='\n'){
			if(rwmb<pdb_a[rwmlt.a].b[rwmlt.b].c[rwma]-1){
				for(ishortb=0;ishortb<pdb_a[rwmlt.a].b[rwmlt.b].c[rwma];ishortb++)
					{ ln[rwmln-1].f[rwma].c[ishortb]='\0'; }
				}
			if(rwma<pdb_a[rwmlt.a].b[rwmlt.b].f-1){
				for(ishorta=rwma;ishorta<pdb_a[rwmlt.a].b[rwmlt.b].f;ishorta++){
				for(ishortb=0;ishortb<pdb_a[rwmlt.a].b[rwmlt.b].c[ishorta];ishortb++){
					ln[rwmln-1].f[ishorta].c[ishortb]='\0';
					}}
				}
                        rwmb=pdb_a[rwmlt.a].b[rwmlt.b].c[rwma];
                        rwma=pdb_a[rwmlt.a].b[rwmlt.b].f;
                        }
                if(ln[rwmln-1].f[rwma].c[rwmb]==EOF){
                        read_neek("EOF",(rwmln),(rwma+1));
                        }
		itr++;
                }
        ln[rwmln-1].f[rwma].c[rwmb]='\0';
        } 
//if(DEBUG>=0){printf("rwm_line 13\n");}
/*if((WARN==0)&&(SILENT!=0)){fflush(OUTC);} */
return;
}

