/************************** pdb_def() ******************************/
/* This is part of the program's initialization.  It assigns rule values to
 * the data structures for line formats.  In other words, it tells the program
 * what class each line type is, how many fields the line has, how long each
 * field is and what type of data exists in that field.
 *
 * Data Format :  
 *
 * Record name  : pdb_a[ta].b[tb].typ 	=	 "Name", e.g., "ATOM"
 *
 * Number of fields : pdb_a[ta].b[tb].f 	=	 # of fields 
 *
 * 	NOTE:  This program differs from the pdb manual in that each space
 * 	in the 80-character width is given a field designation.  In the
 * 	manual, some space numbers are assigned to a field.  They are here.
 *
 * Number of characters each field : 
 * 	pdb_a[ta].b[tb].c[Fn] 	=	 # of characters in field Fn 
 *
 * Data type each field : 
 * 	pdb_a[ta].b[tb].t[Fn] 	=	 data type for of characters in field Fn 
 *	
 *	N/A in any comment means this field is not assigned in the pdb rule
 *	set and should be entered as empty spaces in the file.
 *
 *  IN OUTPUT
 *	line numbers and field numbers start counting at 1 (for human ease)
 *	SO....  add '1' to field numbers as listed below to translate
 *
 * Author: Lachele Foley
*/
#include <load_pdb.h>
//#include "../inc/load_pdb.h"
void pdb_def(){
/* Class for this set of Records : Once Only -- Class 0 
 * 	Class 0 records should occur, at most, once in any pdb file */
/* RECORD Name : CRYST1 */
strcpy(pdb_a[0].b[0].typ,"CRYST1");  /*  */
pdb_a[0].b[0].f		=	11;  /*  */
pdb_a[0].b[0].c[0]	=	6;  /* RECORD NAME */
pdb_a[0].b[0].c[1]	=	9;  /* a */
pdb_a[0].b[0].c[2]	=	9;  /* b */
pdb_a[0].b[0].c[3]	=	9;  /* c */
pdb_a[0].b[0].c[4]	=	7;  /* alpha */
pdb_a[0].b[0].c[5]	=	7;  /* beta */
pdb_a[0].b[0].c[6]	=	7;  /* gamma */
pdb_a[0].b[0].c[7]	=	1;  /* N/A */
pdb_a[0].b[0].c[8]	=	11;  /* sGroup */
pdb_a[0].b[0].c[9]	=	4;  /* z */
pdb_a[0].b[0].c[10]	=	10;  /* N/A */
pdb_a[0].b[0].t[0]	=	's';  /* RECORD NAME */
pdb_a[0].b[0].t[1]	=	'f';  /* a */
pdb_a[0].b[0].t[2]	=	'f';  /* b */
pdb_a[0].b[0].t[3]	=	'f';  /* c */
pdb_a[0].b[0].t[4]	=	'f';  /* alpha */
pdb_a[0].b[0].t[5]	=	'f';  /* beta */
pdb_a[0].b[0].t[6]	=	'f';  /* gamma */
pdb_a[0].b[0].t[7]	=	's';  /* N/A */
pdb_a[0].b[0].t[8]	=	's';  /* sGroup */
pdb_a[0].b[0].t[9]	=	'i';  /* z */
pdb_a[0].b[0].t[10]	=	's';  /* N/A */
/* RECORD Name : END */
strcpy(pdb_a[0].b[1].typ,"END");  /*  */
pdb_a[0].b[1].f		=	2;  /*  */
pdb_a[0].b[1].c[0]	=	6;  /* RECORD NAME */
pdb_a[0].b[1].c[1]	=	74;  /* N/A */
pdb_a[0].b[1].t[0]	=	's';  /* RECORD NAME */
pdb_a[0].b[1].t[1]	=	's';  /* N/A */
/* RECORD Name : HEADER */
strcpy(pdb_a[0].b[2].typ,"HEADER");  /*  */
pdb_a[0].b[2].f		=	7;  /*  */
pdb_a[0].b[2].c[0]	=	6;  /* RECORD NAME */
pdb_a[0].b[2].c[1]	=	4;  /* N/A */
pdb_a[0].b[2].c[2]	=	40;  /* classification   */
pdb_a[0].b[2].c[3]	=	9;  /* depDate          */
pdb_a[0].b[2].c[4]	=	3;  /* N/A */
pdb_a[0].b[2].c[5]	=	4;  /* idCode           */
pdb_a[0].b[2].c[6]	=	14;  /* N/A */
pdb_a[0].b[2].t[0]	=	's';  /* RECORD NAME */
pdb_a[0].b[2].t[1]	=	's';  /* N/A */
pdb_a[0].b[2].t[2]	=	's';  /* classification   */
pdb_a[0].b[2].t[3]	=	's';  /* depDate          */
pdb_a[0].b[2].t[4]	=	's';  /* N/A */
pdb_a[0].b[2].t[5]	=	's';  /* idCode           */
pdb_a[0].b[2].t[6]	=	's';  /* N/A */
/* RECORD Name : MASTER */
strcpy(pdb_a[0].b[3].typ,"MASTER");  /*  */
pdb_a[0].b[3].f		=	15;  /*  */
pdb_a[0].b[3].c[0]	=	6;  /* RECORD NAME */
pdb_a[0].b[3].c[1]	=	4;  /* N/A */
pdb_a[0].b[3].c[2]	=	5;  /* numRemark      */
pdb_a[0].b[3].c[3]	=	5;  /* 0 */
pdb_a[0].b[3].c[4]	=	5;  /* numHet         */
pdb_a[0].b[3].c[5]	=	5;  /* numHelix       */
pdb_a[0].b[3].c[6]	=	5;  /* numSheet       */
pdb_a[0].b[3].c[7]	=	5;  /* numTurn        */
pdb_a[0].b[3].c[8]	=	5;  /* numSite        */
pdb_a[0].b[3].c[9]	=	5;  /* numXform       */
pdb_a[0].b[3].c[10]	=	5;  /* numCoord       */
pdb_a[0].b[3].c[11]	=	5;  /* numTer         */
pdb_a[0].b[3].c[12]	=	5;  /* numConect      */
pdb_a[0].b[3].c[13]	=	5;  /* numSeq         */
pdb_a[0].b[3].c[14]	=	10;  /* N/A */
pdb_a[0].b[3].t[0]	=	's';  /* RECORD NAME */
pdb_a[0].b[3].t[1]	=	's';  /* N/A */
pdb_a[0].b[3].t[2]	=	'i';  /* numRemark      */
pdb_a[0].b[3].t[3]	=	'i';  /* 0 */
pdb_a[0].b[3].t[4]	=	'i';  /* numHet         */
pdb_a[0].b[3].t[5]	=	'i';  /* numHelix       */
pdb_a[0].b[3].t[6]	=	'i';  /* numSheet       */
pdb_a[0].b[3].t[7]	=	'i';  /* numTurn        */
pdb_a[0].b[3].t[8]	=	'i';  /* numSite        */
pdb_a[0].b[3].t[9]	=	'i';  /* numXform       */
pdb_a[0].b[3].t[10]	=	'i';  /* numCoord       */
pdb_a[0].b[3].t[11]	=	'i';  /* numTer         */
pdb_a[0].b[3].t[12]	=	'i';  /* numConect      */
pdb_a[0].b[3].t[13]	=	'i';  /* numSeq         */
pdb_a[0].b[3].t[14]	=	's';  /* N/A */
/* RECORD Name : ORIGXn */
strcpy(pdb_a[0].b[4].typ,"ORIGXn");  /*  */
pdb_a[0].b[4].f		=	2;  /*  */
pdb_a[0].b[4].c[0]	=	6;  /* RECORD NAME */
pdb_a[0].b[4].c[1]	=	74;  /* RECORD Text */
pdb_a[0].b[4].t[0]	=	's';  /* RECORD NAME */
pdb_a[0].b[4].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SCALEn */
strcpy(pdb_a[0].b[5].typ,"SCALEn");  /*  */
pdb_a[0].b[5].f		=	2;  /*  */
pdb_a[0].b[5].c[0]	=	6;  /* RECORD NAME */
pdb_a[0].b[5].c[1]	=	74;  /* RECORD Text */
pdb_a[0].b[5].t[0]	=	's';  /* RECORD NAME */
pdb_a[0].b[5].t[1]	=	's';  /* RECORD NAME */
/* Class 1:  Once - many line :  */
/* RECORD Name : AUTHOR */
strcpy(pdb_a[1].b[0].typ,"AUTHOR");  /*  */
pdb_a[1].b[0].f		=	2;  /*  */
pdb_a[1].b[0].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[0].c[1]	=	74;  /* RECORD Text */
pdb_a[1].b[0].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[0].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : CAVEAT */
strcpy(pdb_a[1].b[1].typ,"CAVEAT");  /*  */
pdb_a[1].b[1].f		=	2;  /*  */
pdb_a[1].b[1].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[1].c[1]	=	74;  /* RECORD Text */
pdb_a[1].b[1].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[1].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : COMPND */
strcpy(pdb_a[1].b[2].typ,"COMPND");  /*  */
pdb_a[1].b[2].f		=	2;  /*  */
pdb_a[1].b[2].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[2].c[1]	=	74;  /* RECORD Text */
pdb_a[1].b[2].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[2].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : EXPDTA */
strcpy(pdb_a[1].b[3].typ,"EXPDTA");  /*  */
pdb_a[1].b[3].f		=	2;  /*  */
pdb_a[1].b[3].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[3].c[1]	=	74;  /* RECORD Text */
pdb_a[1].b[3].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[3].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : KEYWDS */
strcpy(pdb_a[1].b[4].typ,"KEYWDS");  /*  */
pdb_a[1].b[4].f		=	2;  /*  */
pdb_a[1].b[4].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[4].c[1]	=	74;  /* RECORD Text */
pdb_a[1].b[4].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[4].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : OBSLTE */
strcpy(pdb_a[1].b[5].typ,"OBSLTE");  /*  */
pdb_a[1].b[5].f		=	2;  /*  */
pdb_a[1].b[5].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[5].c[1]	=	74;  /* RECORD Text */
pdb_a[1].b[5].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[5].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SOURCE */
strcpy(pdb_a[1].b[6].typ,"SOURCE");  /*  */
pdb_a[1].b[6].f		=	2;  /*  */
pdb_a[1].b[6].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[6].c[1]	=	74;  /* RECORD Text */
pdb_a[1].b[6].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[6].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SPRSDE */
strcpy(pdb_a[1].b[7].typ,"SPRSDE");  /*  */
pdb_a[1].b[7].f		=	2;  /*  */
pdb_a[1].b[7].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[7].c[1]	=	74;  /* RECORD Text */
pdb_a[1].b[7].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[7].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : TITLE */
strcpy(pdb_a[1].b[8].typ,"TITLE");  /*  */
pdb_a[1].b[8].f		=	5;  /*  */
pdb_a[1].b[8].c[0]	=	6;  /* RECORD NAME */
pdb_a[1].b[8].c[1]	=	2;  /* N/A */
pdb_a[1].b[8].c[2]	=	2;  /* continuation    */
pdb_a[1].b[8].c[3]	=	60;  /* title           */
pdb_a[1].b[8].c[4]	=	10;  /* N/A */
pdb_a[1].b[8].t[0]	=	's';  /* RECORD NAME */
pdb_a[1].b[8].t[1]	=	's';  /* N/A */
pdb_a[1].b[8].t[2]	=	's';  /* continuation    */
pdb_a[1].b[8].t[3]	=	's';  /* title           */
pdb_a[1].b[8].t[4]	=	's';  /* N/A */
/* Class 2:  Many - One line each :  */
/* RECORD Name : ANISOU */
strcpy(pdb_a[2].b[0].typ,"ANISOU");  /*  */
pdb_a[2].b[0].f		=	2;  /*  */
pdb_a[2].b[0].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[0].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[0].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[0].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : ATOM */
strcpy(pdb_a[2].b[1].typ,"ATOM");  /*  */
pdb_a[2].b[1].f		=	20;  /*  */
pdb_a[2].b[1].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[1].c[1]	=	5;  /* serial         */
pdb_a[2].b[1].c[2]	=	1;  /* N/A */
pdb_a[2].b[1].c[3]	=	4;  /* name           */
pdb_a[2].b[1].c[4]	=	1;  /* altLoc         */
pdb_a[2].b[1].c[5]	=	3;  /* resName        */
pdb_a[2].b[1].c[6]	=	1;  /* N/A */
pdb_a[2].b[1].c[7]	=	1;  /* chainID        */
pdb_a[2].b[1].c[8]	=	4;  /* resSeq         */
pdb_a[2].b[1].c[9]	=	1;  /* iCode          */
pdb_a[2].b[1].c[10]	=	3;  /* N/A */
pdb_a[2].b[1].c[11]	=	8;  /* x */
pdb_a[2].b[1].c[12]	=	8;  /* y */
pdb_a[2].b[1].c[13]	=	8;  /* z */
pdb_a[2].b[1].c[14]	=	6;  /* occupancy      */
pdb_a[2].b[1].c[15]	=	6;  /* tempFactor     */
pdb_a[2].b[1].c[16]	=	6;  /* N/A */
pdb_a[2].b[1].c[17]	=	4;  /* segID          */
pdb_a[2].b[1].c[18]	=	2;  /* element        */
pdb_a[2].b[1].c[19]	=	2;  /* charge         */
pdb_a[2].b[1].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[1].t[1]	=	'i';  /* serial         */
pdb_a[2].b[1].t[2]	=	's';  /* N/A */
pdb_a[2].b[1].t[3]	=	's';  /* name           */
pdb_a[2].b[1].t[4]	=	'c';  /* altLoc         */
pdb_a[2].b[1].t[5]	=	's';  /* resName        */
pdb_a[2].b[1].t[6]	=	's';  /* N/A */
pdb_a[2].b[1].t[7]	=	'c';  /* chainID        */
pdb_a[2].b[1].t[8]	=	'i';  /* resSeq         */
pdb_a[2].b[1].t[9]	=	'c';  /* iCode          */
pdb_a[2].b[1].t[10]	=	's';  /* N/A */
pdb_a[2].b[1].t[11]	=	'f';  /* x */
pdb_a[2].b[1].t[12]	=	'f';  /* y */
pdb_a[2].b[1].t[13]	=	'f';  /* z */
pdb_a[2].b[1].t[14]	=	'f';  /* occupancy      */
pdb_a[2].b[1].t[15]	=	'f';  /* tempFactor     */
pdb_a[2].b[1].t[16]	=	's';  /* N/A */
pdb_a[2].b[1].t[17]	=	's';  /* segID          */
pdb_a[2].b[1].t[18]	=	's';  /* element        */
pdb_a[2].b[1].t[19]	=	's';  /* charge         */
/* RECORD Name : CISPEP */
strcpy(pdb_a[2].b[2].typ,"CISPEP");  /*  */
pdb_a[2].b[2].f		=	2;  /*  */
pdb_a[2].b[2].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[2].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[2].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[2].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : CONECT */
strcpy(pdb_a[2].b[3].typ,"CONECT");  /*  */
pdb_a[2].b[3].f		=	2;  /*  */
pdb_a[2].b[3].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[3].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[3].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[3].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : DBREF */
strcpy(pdb_a[2].b[4].typ,"DBREF");  /*  */
pdb_a[2].b[4].f		=	2;  /*  */
pdb_a[2].b[4].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[4].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[4].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[4].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : HELIX */
strcpy(pdb_a[2].b[5].typ,"HELIX");  /*  */
pdb_a[2].b[5].f		=	2;  /*  */
pdb_a[2].b[5].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[5].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[5].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[5].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : HET */
strcpy(pdb_a[2].b[6].typ,"HET");  /*  */
pdb_a[2].b[6].f		=	12;  /*  */
pdb_a[2].b[6].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[6].c[1]	=	1;  /* N/A */
pdb_a[2].b[6].c[2]	=	3;  /* hetID          */
pdb_a[2].b[6].c[3]	=	2;  /* N/A */
pdb_a[2].b[6].c[4]	=	1;  /* ChainID        */
pdb_a[2].b[6].c[5]	=	4;  /* seqNum         */
pdb_a[2].b[6].c[6]	=	1;  /* iCode          */
pdb_a[2].b[6].c[7]	=	2;  /* N/A */
pdb_a[2].b[6].c[8]	=	5;  /* numHetAtoms    */
pdb_a[2].b[6].c[9]	=	5;  /* N/A */
pdb_a[2].b[6].c[10]	=	40;  /* text           */
pdb_a[2].b[6].c[11]	=	10;  /* N/A */
pdb_a[2].b[6].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[6].t[1]	=	's';  /* N/A */
pdb_a[2].b[6].t[2]	=	's';  /* hetID          */
pdb_a[2].b[6].t[3]	=	's';  /* N/A */
pdb_a[2].b[6].t[4]	=	'c';  /* ChainID        */
pdb_a[2].b[6].t[5]	=	'i';  /* seqNum         */
pdb_a[2].b[6].t[6]	=	'c';  /* iCode          */
pdb_a[2].b[6].t[7]	=	's';  /* N/A */
pdb_a[2].b[6].t[8]	=	'i';  /* numHetAtoms    */
pdb_a[2].b[6].t[9]	=	's';  /* N/A */
pdb_a[2].b[6].t[10]	=	's';  /* text           */
pdb_a[2].b[6].t[11]	=	's';  /* N/A */
/* RECORD Name : HETSYN */
strcpy(pdb_a[2].b[7].typ,"HETSYN");  /*  */
pdb_a[2].b[7].f		=	8;  /*  */
pdb_a[2].b[7].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[7].c[1]	=	2;  /* N/A */
pdb_a[2].b[7].c[2]	=	2;  /* continuation    */
pdb_a[2].b[7].c[3]	=	1;  /* N/A */
pdb_a[2].b[7].c[4]	=	3;  /* hetID           */
pdb_a[2].b[7].c[5]	=	1;  /* N/A */
pdb_a[2].b[7].c[6]	=	55;  /* hetSynonyms     */
pdb_a[2].b[7].c[7]	=	10;  /* N/A */
pdb_a[2].b[7].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[7].t[1]	=	's';  /* N/A */
pdb_a[2].b[7].t[2]	=	's';  /* continuation    */
pdb_a[2].b[7].t[3]	=	's';  /* N/A */
pdb_a[2].b[7].t[4]	=	's';  /* hetID           */
pdb_a[2].b[7].t[5]	=	's';  /* N/A */
pdb_a[2].b[7].t[6]	=	's';  /* hetSynonyms     */
pdb_a[2].b[7].t[7]	=	's';  /* N/A */
/* RECORD Name : HYDBND */
strcpy(pdb_a[2].b[8].typ,"HYDBND");  /*  */
pdb_a[2].b[8].f		=	2;  /*  */
pdb_a[2].b[8].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[8].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[8].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[8].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : LINK */
strcpy(pdb_a[2].b[9].typ,"LINK");  /*  */
pdb_a[2].b[9].f		=	2;  /*  */
pdb_a[2].b[9].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[9].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[9].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[9].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : MODRES */
strcpy(pdb_a[2].b[10].typ,"MODRES");  /*  */
pdb_a[2].b[10].f		=	2;  /*  */
pdb_a[2].b[10].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[10].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[10].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[10].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : MTRIXn */
strcpy(pdb_a[2].b[11].typ,"MTRIXn");  /*  */
pdb_a[2].b[11].f		=	2;  /*  */
pdb_a[2].b[11].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[11].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[11].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[11].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : REVDAT */
strcpy(pdb_a[2].b[12].typ,"REVDAT");  /*  */
pdb_a[2].b[12].f		=	2;  /*  */
pdb_a[2].b[12].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[12].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[12].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[12].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SEQADV */
strcpy(pdb_a[2].b[13].typ,"SEQADV");  /*  */
pdb_a[2].b[13].f		=	2;  /*  */
pdb_a[2].b[13].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[13].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[13].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[13].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SEQRES */
strcpy(pdb_a[2].b[14].typ,"SEQRES");  /*  */
pdb_a[2].b[14].f		=	2;  /*  */
pdb_a[2].b[14].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[14].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[14].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[14].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SHEET */
strcpy(pdb_a[2].b[15].typ,"SHEET");  /*  */
pdb_a[2].b[15].f		=	2;  /*  */
pdb_a[2].b[15].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[15].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[15].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[15].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SIGATM */
strcpy(pdb_a[2].b[16].typ,"SIGATM");  /*  */
pdb_a[2].b[16].f		=	2;  /*  */
pdb_a[2].b[16].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[16].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[16].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[16].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SIGUIJ */
strcpy(pdb_a[2].b[17].typ,"SIGUIJ");  /*  */
pdb_a[2].b[17].f		=	2;  /*  */
pdb_a[2].b[17].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[17].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[17].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[17].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SITE */
strcpy(pdb_a[2].b[18].typ,"SITE");  /*  */
pdb_a[2].b[18].f		=	2;  /*  */
pdb_a[2].b[18].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[18].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[18].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[18].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SLTBRG */
strcpy(pdb_a[2].b[19].typ,"SLTBRG");  /*  */
pdb_a[2].b[19].f		=	2;  /*  */
pdb_a[2].b[19].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[19].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[19].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[19].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : SSBOND */
strcpy(pdb_a[2].b[20].typ,"SSBOND");  /*  */
pdb_a[2].b[20].f		=	2;  /*  */
pdb_a[2].b[20].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[20].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[20].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[20].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : TURN */
strcpy(pdb_a[2].b[21].typ,"TURN");  /*  */
pdb_a[2].b[21].f		=	2;  /*  */
pdb_a[2].b[21].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[21].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[21].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[21].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : TVECT */
strcpy(pdb_a[2].b[22].typ,"TVECT");  /*  */
pdb_a[2].b[22].f		=	2;  /*  */
pdb_a[2].b[22].c[0]	=	6;  /* RECORD NAME */
pdb_a[2].b[22].c[1]	=	74;  /* RECORD Text */
pdb_a[2].b[22].t[0]	=	's';  /* RECORD NAME */
pdb_a[2].b[22].t[1]	=	's';  /* RECORD NAME */
/* Class 3:  Many -- many line each :  */
/* RECORD Name : FORMUL */
strcpy(pdb_a[3].b[0].typ,"FORMUL");  /*  */
pdb_a[3].b[0].f		=	2;  /*  */
pdb_a[3].b[0].c[0]	=	6;  /* RECORD NAME */
pdb_a[3].b[0].c[1]	=	74;  /* RECORD Text */
pdb_a[3].b[0].t[0]	=	's';  /* RECORD NAME */
pdb_a[3].b[0].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : HETATM */
strcpy(pdb_a[3].b[1].typ,"HETATM");  /*  */
pdb_a[3].b[1].f		=	20;  /*  */
pdb_a[3].b[1].c[0]	=	6;  /* RECORD NAME */
pdb_a[3].b[1].c[1]	=	5;  /* serial         */
pdb_a[3].b[1].c[2]	=	1;  /* N/A */
pdb_a[3].b[1].c[3]	=	4;  /* name           */
pdb_a[3].b[1].c[4]	=	1;  /* altLoc         */
pdb_a[3].b[1].c[5]	=	3;  /* resName        */
pdb_a[3].b[1].c[6]	=	1;  /* N/A */
pdb_a[3].b[1].c[7]	=	1;  /* chainID        */
pdb_a[3].b[1].c[8]	=	4;  /* resSeq         */
pdb_a[3].b[1].c[9]	=	1;  /* iCode          */
pdb_a[3].b[1].c[10]	=	3;  /* N/A */
pdb_a[3].b[1].c[11]	=	8;  /* x */
pdb_a[3].b[1].c[12]	=	8;  /* y */
pdb_a[3].b[1].c[13]	=	8;  /* z */
pdb_a[3].b[1].c[14]	=	6;  /* occupancy      */
pdb_a[3].b[1].c[15]	=	6;  /* tempFactor     */
pdb_a[3].b[1].c[16]	=	6;  /* N/A */
pdb_a[3].b[1].c[17]	=	4;  /* segID          */
pdb_a[3].b[1].c[18]	=	2;  /* element        */
pdb_a[3].b[1].c[19]	=	2;  /* charge         */
pdb_a[3].b[1].t[0]	=	's';  /* RECORD NAME */
pdb_a[3].b[1].t[1]	=	'i';  /* serial         */
pdb_a[3].b[1].t[2]	=	's';  /* N/A */
pdb_a[3].b[1].t[3]	=	's';  /* name           */
pdb_a[3].b[1].t[4]	=	'c';  /* altLoc         */
pdb_a[3].b[1].t[5]	=	's';  /* resName        */
pdb_a[3].b[1].t[6]	=	's';  /* N/A */
pdb_a[3].b[1].t[7]	=	'c';  /* chainID        */
pdb_a[3].b[1].t[8]	=	'i';  /* resSeq         */
pdb_a[3].b[1].t[9]	=	'c';  /* iCode          */
pdb_a[3].b[1].t[10]	=	's';  /* N/A */
pdb_a[3].b[1].t[11]	=	'f';  /* x */
pdb_a[3].b[1].t[12]	=	'f';  /* y */
pdb_a[3].b[1].t[13]	=	'f';  /* z */
pdb_a[3].b[1].t[14]	=	'f';  /* occupancy      */
pdb_a[3].b[1].t[15]	=	'f';  /* tempFactor     */
pdb_a[3].b[1].t[16]	=	's';  /* N/A */
pdb_a[3].b[1].t[17]	=	's';  /* segID          */
pdb_a[3].b[1].t[18]	=	's';  /* element        */
pdb_a[3].b[1].t[19]	=	's';  /* charge         */
/* RECORD Name : HETNAM */
strcpy(pdb_a[3].b[2].typ,"HETNAM");  /*  */
pdb_a[3].b[2].f		=	8;  /*  */
pdb_a[3].b[2].c[0]	=	6;  /* RECORD NAME */
pdb_a[3].b[2].c[1]	=	2;  /* N/A */
pdb_a[3].b[2].c[2]	=	2;  /* continuation    */
pdb_a[3].b[2].c[3]	=	1;  /* N/A */
pdb_a[3].b[2].c[4]	=	3;  /* hetID           */
pdb_a[3].b[2].c[5]	=	1;  /* N/A */
pdb_a[3].b[2].c[6]	=	55;  /* text            */
pdb_a[3].b[2].c[7]	=	10;  /* N/A */
pdb_a[3].b[2].t[0]	=	's';  /* RECORD NAME */
pdb_a[3].b[2].t[1]	=	's';  /* N/A */
pdb_a[3].b[2].t[2]	=	's';  /* continuation    */
pdb_a[3].b[2].t[3]	=	's';  /* N/A */
pdb_a[3].b[2].t[4]	=	's';  /* hetID           */
pdb_a[3].b[2].t[5]	=	's';  /* N/A */
pdb_a[3].b[2].t[6]	=	's';  /* text            */
pdb_a[3].b[2].t[7]	=	's';  /* N/A */
/* grouping lines :  */
/* RECORD Name : ENDMDL */
strcpy(pdb_a[4].b[0].typ,"ENDMDL");  /*  */
pdb_a[4].b[0].f		=	2;  /*  */
pdb_a[4].b[0].c[0]	=	6;  /* RECORD NAME */
pdb_a[4].b[0].c[1]	=	74;  /* N/A */
pdb_a[4].b[0].t[0]	=	's';  /* RECORD NAME */
pdb_a[4].b[0].t[1]	=	's';  /* N/A */
/* RECORD Name : MODEL */
strcpy(pdb_a[4].b[1].typ,"MODEL");  /*  */
pdb_a[4].b[1].f		=	4;  /*  */
pdb_a[4].b[1].c[0]	=	6;  /* RECORD NAME */
pdb_a[4].b[1].c[1]	=	4;  /* N/A */
pdb_a[4].b[1].c[2]	=	4;  /* serial         */
pdb_a[4].b[1].c[3]	=	66;  /* N/A */
pdb_a[4].b[1].t[0]	=	's';  /* RECORD NAME */
pdb_a[4].b[1].t[1]	=	's';  /* N/A */
pdb_a[4].b[1].t[2]	=	'i';  /* serial         */
pdb_a[4].b[1].t[3]	=	's';  /* N/A */
/* RECORD Name : TER */
strcpy(pdb_a[4].b[2].typ,"TER");  /*  */
pdb_a[4].b[2].f		=	9;  /*  */
pdb_a[4].b[2].c[0]	=	6;  /* RECORD NAME */
pdb_a[4].b[2].c[1]	=	5;  /* serial        */
pdb_a[4].b[2].c[2]	=	6;  /* N/A */
pdb_a[4].b[2].c[3]	=	3;  /* resName       */
pdb_a[4].b[2].c[4]	=	1;  /* N/A */
pdb_a[4].b[2].c[5]	=	1;  /* chainID       */
pdb_a[4].b[2].c[6]	=	4;  /* resSeq        */
pdb_a[4].b[2].c[7]	=	1;  /* iCode         */
pdb_a[4].b[2].c[8]	=	53;  /* N/A */
pdb_a[4].b[2].t[0]	=	's';  /* RECORD NAME */
pdb_a[4].b[2].t[1]	=	'i';  /* serial        */
pdb_a[4].b[2].t[2]	=	's';  /* N/A */
pdb_a[4].b[2].t[3]	=	's';  /* resName       */
pdb_a[4].b[2].t[4]	=	's';  /* N/A */
pdb_a[4].b[2].t[5]	=	'c';  /* chainID       */
pdb_a[4].b[2].t[6]	=	'i';  /* resSeq        */
pdb_a[4].b[2].t[7]	=	'c';  /* iCode         */
pdb_a[4].b[2].t[8]	=	's';  /* N/A */
/* Special :  */
/* RECORD Name : JRNL */
strcpy(pdb_a[5].b[0].typ,"JRNL");  /*  */
pdb_a[5].b[0].f		=	2;  /*  */
pdb_a[5].b[0].c[0]	=	6;  /* RECORD NAME */
pdb_a[5].b[0].c[1]	=	74;  /* RECORD Text */
pdb_a[5].b[0].t[0]	=	's';  /* RECORD NAME */
pdb_a[5].b[0].t[1]	=	's';  /* RECORD NAME */
/* RECORD Name : REMARK */
strcpy(pdb_a[5].b[1].typ,"REMARK");  /*  */
pdb_a[5].b[1].f		=	6;  /*  */
pdb_a[5].b[1].c[0]	=	6;  /* RECORD NAME */
pdb_a[5].b[1].c[1]	=	1;  /* N/A */
pdb_a[5].b[1].c[2]	=	3;  /* remarkNum       */
pdb_a[5].b[1].c[3]	=	1;  /* N/A */
pdb_a[5].b[1].c[4]	=	59;  /* empty           */
pdb_a[5].b[1].c[5]	=	10;  /* N/A */
pdb_a[5].b[1].t[0]	=	's';  /* RECORD NAME */
pdb_a[5].b[1].t[1]	=	's';  /* N/A */
pdb_a[5].b[1].t[2]	=	'i';  /* remarkNum       */
pdb_a[5].b[1].t[3]	=	's';  /* N/A */
pdb_a[5].b[1].t[4]	=	's';  /* empty           */
pdb_a[5].b[1].t[5]	=	's';  /* N/A */
return;
}
