/**************************  get_type() *******************************/
/* get_type(mca) finds the entry type for that line based on the first word */
/* Author: Lachele Foley*/
#include <load_pdb.h>
//#include "../inc/load_pdb.h"
linetype get_type(char* gca){ 
int gta=0;
char gt_other[20];
linetype gtlt;
if(DEBUG>=1){printf("gca is >>>%s<<<\n", gca);}
gtlt.a=-1;
gtlt.b=-1;

/* CLASS 0 */
if(strcmp(gca,"CRYST1")==0){
	gtlt.a=0;
	gtlt.b=0;
	}  /* CRYST1 */
if(strcmp(gca,"END   ")==0){
	gtlt.a=0;
	gtlt.b=1;
	}  /* END */
if(strcmp(gca,"HEADER")==0){
	gtlt.a=0;
	gtlt.b=2;
	}  /* HEADER */
if(strcmp(gca,"MASTER")==0){
	gtlt.a=0;
	gtlt.b=3;
	}  /* MASTER */
/* In case ORIGXn and SCALEn are taken literally: */
if(strcmp(gca,"ORIGXn")==0){
	gtlt.a=0;
	gtlt.b=4;
	}  /* ORIGXn */
if(strcmp(gca,"SCALEn")==0){
	gtlt.a=0;
	gtlt.b=5;
	}  /* SCALEn */
/* In case ORIGXn and SCALEn are given n=variable: */
for(gta=0;gta<5;gta++){
	gt_other[gta]=gca[gta];
	}
gt_other[5]='\0';
if(strcmp(gt_other,"ORIGX")==0){
	gtlt.a=0;
	gtlt.b=4;
	}  /* ORIGXn */ 
if(strcmp(gt_other,"SCALE")==0){
	gtlt.a=0;
	gtlt.b=5;
	}  /* ORIGXn */ 

/* CLASS 1 */
if(strcmp(gca,"AUTHOR")==0){
	gtlt.a=1;
	gtlt.b=0;
	}  /* AUTHOR */
if(strcmp(gca,"CAVEAT")==0){
	gtlt.a=1;
	gtlt.b=1;
	}  /* CAVEAT */
if(strcmp(gca,"COMPND")==0){
	gtlt.a=1;
	gtlt.b=2;
	}  /* COMPND */
if(strcmp(gca,"EXPDTA")==0){
	gtlt.a=1;
	gtlt.b=3;
	}  /* EXPDTA */
if(strcmp(gca,"KEYWDS")==0){
	gtlt.a=1;
	gtlt.b=4;
	}  /* KEYWDS */
if(strcmp(gca,"OBSLTE")==0){
	gtlt.a=1;
	gtlt.b=5;
	}  /* OBSLTE */
if(strcmp(gca,"SOURCE")==0){
	gtlt.a=1;
	gtlt.b=6;
	}  /* SOURCE */
if(strcmp(gca,"SPRSDE")==0){
	gtlt.a=1;
	gtlt.b=7;
	}  /* SPRSDE */
if(strcmp(gca,"TITLE ")==0){
	gtlt.a=1;
	gtlt.b=8;
	}  /* TITLE */

/* CLASS 2 */
if(strcmp(gca,"ANISOU")==0){
	gtlt.a=2;
	gtlt.b=0;
	}  /* ANISOU */
if(strcmp(gca,"ATOM  ")==0){
	gtlt.a=2;
	gtlt.b=1;
	}  /* ATOM */
if(strcmp(gca,"CISPEP")==0){
	gtlt.a=2;
	gtlt.b=2;
	}  /* CISPEP */
if(strcmp(gca,"CONECT")==0){
	gtlt.a=2;
	gtlt.b=3;
	}  /* CONECT */
if(strcmp(gca,"DBREF ")==0){
	gtlt.a=2;
	gtlt.b=4;
	}  /* DBREF */
if(strcmp(gca,"HELIX ")==0){
	gtlt.a=2;
	gtlt.b=5;
	}  /* HELIX */
if(strcmp(gca,"HET   ")==0){
	gtlt.a=2;
	gtlt.b=6;
	}  /* HET */
if(strcmp(gca,"HETSYN")==0){
	gtlt.a=2;
	gtlt.b=7;
	}  /* HETSYN */
if(strcmp(gca,"HYDBND")==0){
	gtlt.a=2;
	gtlt.b=8;
	}  /* HYDBND */
if(strcmp(gca,"LINK  ")==0){
	gtlt.a=2;
	gtlt.b=9;
	}  /* LINK */
if(strcmp(gca,"MODRES")==0){
	gtlt.a=2;
	gtlt.b=10;
	}  /* MODRES */
if(strcmp(gca,"MTRIXn")==0){
	gtlt.a=2;
	gtlt.b=11;
	}  /* MTRIXn */
if(strcmp(gca,"REVDAT")==0){
	gtlt.a=2;
	gtlt.b=12;
	}  /* REVDAT */
if(strcmp(gca,"SEQADV")==0){
	gtlt.a=2;
	gtlt.b=13;
	}  /* SEQADV */
if(strcmp(gca,"SEQRES")==0){
	gtlt.a=2;
	gtlt.b=14;
	}  /* SEQRES */
if(strcmp(gca,"SHEET ")==0){
	gtlt.a=2;
	gtlt.b=15;
	}  /* SHEET */
if(strcmp(gca,"SIGATM")==0){
	gtlt.a=2;
	gtlt.b=16;
	}  /* SIGATM */
if(strcmp(gca,"SIGUIJ")==0){
	gtlt.a=2;
	gtlt.b=17;
	}  /* SIGUIJ */
if(strcmp(gca,"SITE  ")==0){
	gtlt.a=2;
	gtlt.b=18;
	}  /* SITE */
if(strcmp(gca,"SLTBRG")==0){
	gtlt.a=2;
	gtlt.b=19;
	}  /* SLTBRG */
if(strcmp(gca,"SSBOND")==0){
	gtlt.a=2;
	gtlt.b=20;
	}  /* SSBOND */
if(strcmp(gca,"TURN  ")==0){
	gtlt.a=2;
	gtlt.b=21;
	}  /* TURN */
if(strcmp(gca,"TVECT ")==0){
	gtlt.a=2;
	gtlt.b=22;
	}  /* TVECT */

/* CLASS 3 */
if(strcmp(gca,"FORMUL")==0){
	gtlt.a=3;
	gtlt.b=0;
	}  /* FORMUL */
if(strcmp(gca,"HETATM")==0){
	gtlt.a=3;
	gtlt.b=1;
	}  /* HETATM */
if(strcmp(gca,"HETNAM")==0){
	gtlt.a=3;
	gtlt.b=2;
	}  /* HETNAM */

/* CLASS 4 */
if(strcmp(gca,"ENDMDL")==0){
	gtlt.a=4;
	gtlt.b=0;
	}  /* ENDMDL */
if(strcmp(gca,"MODEL ")==0){
	gtlt.a=4;
	gtlt.b=1;
	}  /* MODEL */
if(strcmp(gca,"TER   ")==0){
	gtlt.a=4;
	gtlt.b=2;
	}  /* TER */

/* CLASS 5 */
if(strcmp(gca,"JRNL  ")==0){
	gtlt.a=5;
	gtlt.b=0;
	}  /* JRNL */
if(strcmp(gca,"REMARK")==0){
	gtlt.a=5;
	gtlt.b=1;
	}  /* REMARK */ 
if(gtlt.a<0){
printf("String type of %s found in file.  This type not supported.\n",gca);
printf("    [feel free to write it in!]\n");
printf("Deleting entry from output file.\n"); 
gtlt.a=-1;
gtlt.b=-1;
}
return gtlt;
} 

