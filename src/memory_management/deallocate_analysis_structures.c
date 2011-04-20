#include "../inc/mylib.h"
#include "../inc/molecules.h"
#include "../inc/modes.h"
#include "../inc/stats.h"

/** \file deallocate_analysis_structures.c 
\addtogroup MEMORY_MANAGEMENT
\brief   Dealloction routines for structures relevant to analysis.

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
/********** structures from modes.h ****************/
void deallocateVibaddr(vibaddr *va){
 // free top
 if(va[0].vi != NULL && va[0].vi != 0x0){free(va[0].vi);}
 return ;
}
void deallocateStretch(stretch *s){
 if(s[0].Desc != NULL && s[0].Desc != 0x0){free(s[0].Desc);}
 return ;
}
void deallocateBend(bend *b){
 if(b[0].Desc != NULL && b[0].Desc != 0x0){free(b[0].Desc);}
 return ;
}
void deallocateTorsion(torsion *t){
 if(t[0].Desc != NULL && t[0].Desc != 0x0){free(t[0].Desc);}
 return ;
}
void deallocateRingmotion(ringmotion *rm){
 // free top
 if(rm[0].Desc != NULL && rm[0].Desc != 0x0){free(rm[0].Desc);}
 if(rm[0].mi != NULL && rm[0].mi != 0x0){free(rm[0].mi);}
 if(rm[0].c != NULL && rm[0].c != 0x0){free(rm[0].c);}
 return ;
}
void deallocateVibmode(vibmode *vm){
 int i;
 // deallocate each
 if(vm[0].s != NULL && vm[0].s != 0x0){
 	for(i=0;i<vm[0].ns;i++){deallocateStretch(&vm[0].s[i]);}
 	free(vm[0].s);}
 if(vm[0].b != NULL && vm[0].b != 0x0){
 	for(i=0;i<vm[0].nb;i++){deallocateBend(&vm[0].b[i]);}
 	free(vm[0].b);}
 if(vm[0].t != NULL && vm[0].t != 0x0){
 	for(i=0;i<vm[0].nt;i++){deallocateTorsion(&vm[0].t[i]);}
 	free(vm[0].t);}
 if(vm[0].r != NULL && vm[0].r != 0x0){
 	for(i=0;i<vm[0].nr;i++){deallocateRingmotion(&vm[0].r[i]);}
 	free(vm[0].r);}
 // free top
 if(vm[0].Desc != NULL && vm[0].Desc != 0x0){free(vm[0].Desc);}
 return ;
}
void deallocateTwoassign(twoassign *ta){
 // free top
 if(ta[0].s2 != NULL && ta[0].s2 != 0x0){free(ta[0].s2);}
 if(ta[0].b2 != NULL && ta[0].b2 != 0x0){free(ta[0].b2);}
 if(ta[0].t2 != NULL && ta[0].t2 != 0x0){free(ta[0].t2);}
 if(ta[0].s2s != NULL && ta[0].s2s != 0x0){free(ta[0].s2s);}
 if(ta[0].b2s != NULL && ta[0].b2s != 0x0){free(ta[0].b2s);}
 if(ta[0].t2s != NULL && ta[0].t2s != 0x0){free(ta[0].t2s);}
 if(ta[0].s2a != NULL && ta[0].s2a != 0x0){free(ta[0].s2a);}
 if(ta[0].b2a != NULL && ta[0].b2a != 0x0){free(ta[0].b2a);}
 if(ta[0].t2a != NULL && ta[0].t2a != 0x0){free(ta[0].t2a);}
 return ;
}
void deallocateRingassign(ringassign *ra){
 // free top
 if(ra[0].opsA != NULL && ra[0].opsA != 0x0){free(ra[0].opsA);}
 if(ra[0].opss != NULL && ra[0].opss != 0x0){free(ra[0].opss);}
 if(ra[0].opsa != NULL && ra[0].opsa != 0x0){free(ra[0].opsa);}
 if(ra[0].opdA != NULL && ra[0].opdA != 0x0){free(ra[0].opdA);}
 if(ra[0].opdA1G != NULL && ra[0].opdA1G != 0x0){free(ra[0].opdA1G);}
 if(ra[0].opdOth != NULL && ra[0].opdOth != 0x0){free(ra[0].opdOth);}
 if(ra[0].rbA != NULL && ra[0].rbA != 0x0){free(ra[0].rbA);}
 if(ra[0].rbs != NULL && ra[0].rbs != 0x0){free(ra[0].rbs);}
 if(ra[0].rba != NULL && ra[0].rba != 0x0){free(ra[0].rba);}
 if(ra[0].roA != NULL && ra[0].roA != 0x0){free(ra[0].roA);}
 if(ra[0].ros != NULL && ra[0].ros != 0x0){free(ra[0].ros);}
 if(ra[0].roa != NULL && ra[0].roa != 0x0){free(ra[0].roa);}
 if(ra[0].rpA != NULL && ra[0].rpA != 0x0){free(ra[0].rpA);}
 if(ra[0].rps != NULL && ra[0].rps != 0x0){free(ra[0].rps);}
 if(ra[0].rpa != NULL && ra[0].rpa != 0x0){free(ra[0].rpa);}
 if(ra[0].rso != NULL && ra[0].rso != 0x0){free(ra[0].rso);}
 return ;
}
void deallocateAssignment(assignment *a){
 int i;
 // deallocate each
 if(a[0].r != NULL && a[0].r != 0x0){
 	for(i=0;i<a[0].nr;i++){deallocateRingassign(&a[0].r[i]);}
 	free(a[0].r);}
 if(a[0].two != NULL && a[0].two != 0x0){
 	for(i=0;i<a[0].ntwo;i++){deallocateTwoassign(&a[0].two[i]);}
 	free(a[0].two);}
 // free top
 if(a[0].Desc != NULL && a[0].Desc != 0x0){free(a[0].Desc);}
 if(a[0].s != NULL && a[0].s != 0x0){free(a[0].s);}
 if(a[0].b != NULL && a[0].b != 0x0){free(a[0].b);}
 if(a[0].t != NULL && a[0].t != 0x0){free(a[0].t);}
 if(a[0].w != NULL && a[0].w != 0x0){free(a[0].w);}
 if(a[0].rk != NULL && a[0].rk != 0x0){free(a[0].rk);}
 return ;
}
void deallocateAssnBrief(assn_brief *ab){
 if(ab[0].Desc != NULL && ab[0].Desc != 0x0){free(ab[0].Desc);}
 return ;
}
void deallocateAtommode(atommode *a){
 if(a[0].s != NULL && a[0].s != 0x0){free(a[0].s);}
 if(a[0].b != NULL && a[0].b != 0x0){free(a[0].b);}
 if(a[0].t != NULL && a[0].t != 0x0){free(a[0].t);}
 return ;
}
/********** structures from stats.h ****************/
void deallocateStatsarray(statsarray *s){
 if(s[0].d != NULL && s[0].d != 0x0){free(s[0].d);}
 return ;
}
void deallocateAutocorr(autocorr *a){
 if(a[0].a != NULL && a[0].a != 0x0){free(a[0].a);}
 return ;
} 
