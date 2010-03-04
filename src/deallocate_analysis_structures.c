#include "../inc/mylib.h"
#include "../inc/molecules.h"

/* \file deallocate_amber_structures.c 
\addtogroup MEMORY_MANAGEMENT
\brief   Dealloction routines for structures relevant to AMBER.

Begun on 20100304 by BLFoley 

Notes regarding these functions:

* If a structure does not contain any pointers, it can just be freed.
* Structures that contain pointers must have a deallocation function.
* Single-pointers to arrays of simple types (e.g. int) may be freed.
* Single-pointers to arrays of structures might need individual deallocation.
* Double-pointers
	* Be very careful before freeing double-pointed structures
	* 	-- they might point somewhere you don't want freed
	* Freeing the top-level pointer should not interfere with data below 
*/
/*
void deallocateX(X *a){
 int i;
 //set null
 a[0].t = 0x0;
 // deallocate each
 if(a[0].b != NULL && a[0].b != 0x0){
 	for(i=0;i<a[0].nb;i++){deallocateX(&a[0].b[i]);}
 	free(a[0].b);}
 // free each
 if(a[0].OD != NULL && a[0].OD != 0x0){
 	for(i=0;i<a[0].nOD;i++){if(a[0].OD[i] != NULL && a[0].OD[i] != 0x0){free(a[0].OD[i]);}}
	free(a[0].OD);}
 // free top
 if(a[0].N != NULL && a[0].N != 0x0){free(a[0].N);}
 // check and warn if non-null
 if(a[0].VP!= NULL && a[0].VP != 0x0){fprintf(stderr,"WARNING: deallocating X structure with non-NULL void pointer!\n");}
 return ;
}
*/
