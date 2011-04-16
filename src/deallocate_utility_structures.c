#include "../inc/mylib.h"
#include "../inc/gly_codeutils.h"
#include "../inc/gly_fileutils.h"

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
/********** structures from gly_codeutils.h ****************/
void deallocateLlcharset(llcharset *l){
 if(l[0].T != NULL && l[0].T != 0x0){free(l[0].T);}
 return ;
} 
/********** structures from gly_fileutils.h ****************/
void deallocateGlyKeysvals(gly_keysvals *gk){
 int i;
 // free each
 if(gk[0].K != NULL && gk[0].K != 0x0){
 	for(i=0;i<gk[0].n;i++){if(gk[0].K[i] != NULL && gk[0].K[i] != 0x0){free(gk[0].K[i]);}}
	free(gk[0].K);}
 if(gk[0].V != NULL && gk[0].V != 0x0){
 	for(i=0;i<gk[0].n;i++){if(gk[0].V[i] != NULL && gk[0].V[i] != 0x0){free(gk[0].V[i]);}}
	free(gk[0].V);}
 return ;
}
void deallocateFileset(fileset *f){
 //set null
 f[0].F = 0x0;
 if(f[0].N != NULL && f[0].N != 0x0){free(f[0].N);}
 return ;
}
void deallocateFileslurp(fileslurp *f){
 int i;
 if(f[0].L != NULL && f[0].L != 0x0){
 	for(i=0;i<f[0].n;i++){if(f[0].L[i] != NULL && f[0].L[i] != 0x0){free(f[0].L[i]);}}
	free(f[0].L);}
 return ;
}

