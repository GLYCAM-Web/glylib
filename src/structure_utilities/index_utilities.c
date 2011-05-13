/** \file index_utilities.c
*   Purpose:  Functions associated with internal addresses. 
*   Begin 20110417 BLFoley
* 
*/
#include <mylib.h>
#include <molecules.h>

ensindex copy_moli_to_ensi(molindex moli){
  ensindex ensi;
  ensi.E=ensi.A=-1;
  ensi.i=moli.i;
  ensi.m=moli.m;
  ensi.r=moli.r;
  ensi.a=moli.a;
  return ensi;
  }

char is_consistent_moli_moli(molindex mone, molindex mtwo){
  char answer='y';
  if(mone.i!=mtwo.i) answer='i';
  if(mone.m!=mtwo.m) answer='n';
  if(mone.r!=mtwo.r) answer='n';
  if(mone.a!=mtwo.a) answer='n';
  return answer;
  }

char is_consistent_ensi_ensi(ensindex eone, ensindex etwo){
  char answer='y';
  if(eone.E!=etwo.E) answer='E';
  if(eone.A!=etwo.A) answer='A';
  if(eone.i!=etwo.i) answer='i';
  if(eone.m!=etwo.m) answer='n';
  if(eone.r!=etwo.r) answer='n';
  if(eone.a!=etwo.a) answer='n';
  return answer;
  }


char is_consistent_molbond_molbond(molbond mb1, molbond mb2){
  char answer='y',type='\0';

  answer=is_consistent_moli_moli(mb1.s,mb2.s);
  if(answer=='n') return answer;
  answer=is_consistent_moli_moli(mb1.t,mb2.t);
  if(answer=='n') return answer;

  if(mb1.i!=mb2.i) answer='i';
  if(((mb1.typ==NULL)&&(mb2.typ!=NULL))||((mb1.typ!=NULL)&&(mb2.typ==NULL))) type='t';
  if((mb1.typ!=NULL)&&(mb2.typ!=NULL)) {if((&(mb1.typ[0]))!=(&(mb2.typ[0]))) type='t';}
  if(mb1.o!=mb2.o) type='t';

  if(type=='t')
    {
    if(answer=='i') answer='T';
    if(answer=='y') answer='t';
    }
  
  return answer;
  }
char is_consistent_molbond_molbond_inverse(molbond mb1, molbond mb2){
  char answer='y',type='\0';

  answer=is_consistent_moli_moli(mb1.s,mb2.t);
  if(answer=='n') return answer;
  answer=is_consistent_moli_moli(mb1.t,mb2.s);
  if(answer=='n') return answer;

  if(mb1.i!=mb2.i) answer='i';
  if(((mb1.typ==NULL)&&(mb2.typ!=NULL))||((mb1.typ!=NULL)&&(mb2.typ==NULL))) type='t';
  if((mb1.typ!=NULL)&&(mb2.typ!=NULL)) {if((&(mb1.typ[0]))!=(&(mb2.typ[0]))) type='t';}
  if(mb1.o!=mb2.o) type='t';

  if(type=='t')
    {
    if(answer=='i') answer='T';
    if(answer=='y') answer='t';
    }
  
  return answer;
  }



char is_consistent_ensi_moli(ensindex ensi, molindex moli){
char answer='y';
if(ensi.i!=moli.i) answer='i';
if(ensi.m!=moli.m) answer='n';
if(ensi.r!=moli.r) answer='n';
if(ensi.a!=moli.a) answer='n';
return answer;
}
void set_residue_molindexes(residue *r, int mi, int ri){
int i;
r[0].moli.i=-1;
r[0].moli.m=mi;
r[0].moli.r=ri;
r[0].moli.a=-1;
for(i=0;i<r[0].na;i++)
  {
  r[0].a[i].moli.i=-1;
  r[0].a[i].moli.m=mi;
  r[0].a[i].moli.r=ri;
  r[0].a[i].moli.a=i;
  }
}
void set_molecule_molindexes(molecule *m, int mi){
int i;
m[0].mi=mi;
for(i=0;i<m[0].nr;i++)
  {
  set_residue_molindexes(&m[0].r[i],mi,i);
  }
}
void set_assembly_molindexes(assembly *A){
int i;
for(i=0;i<A[0].nm;i++)
  {
  set_molecule_molindexes(&A[0].m[i][0],i);
  }
}
void set_ensemble_molindexes(ensemble *E){
int i;
for(i=0;i<E[0].nm;i++)
  {
  set_molecule_molindexes(&E[0].m[i],i);
  }
}


