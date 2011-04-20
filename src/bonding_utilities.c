/** \file bonding_utilities.c
*   Purpose:  Functions associated with bonding. 
*   Begin 20110417 BLFoley
*
*   Return values for follow functions:
*        0 : Made it all the way to the end.  All ok.
*       -1 : Made it to an appropriate end in the middle.  Don't continue.
* 
*/
#include <mylib.h>
#include <molecules.h>

/* Move to new file
*/
ensindex copy_moli_to_ensi(molindex moli){
  ensindex ensi;
  ensi.E=ensi.A=-1;
  ensi.i=moli.i;
  ensi.m=moli.m;
  ensi.r=moli.r;
  ensi.a=moli.a;
  return ensi;
  }

/* Move to new file
*/
char is_consistent_moli_moli(molindex mone, molindex mtwo){
  char answer='y';
  if(mone.i!=mtwo.i) answer='i';
  if(mone.m!=mtwo.m) answer='n';
  if(mone.r!=mtwo.r) answer='n';
  if(mone.a!=mtwo.a) answer='n';
  return answer;
  }

/* Move to new file
*/
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
/* Move to new file
*/
char is_consistent_molbond_molbond(molbond mb1, molbond mb2){
  char answer='y',type='\0';

  answer=is_consistent_moli_moli(mb1.s,mb2.s);
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

/* Move to new file
*/
char is_consistent_ensi_moli(ensindex ensi, molindex moli){
char answer='y';
if(ensi.i!=moli.i) answer='i';
if(ensi.m!=moli.m) answer='n';
if(ensi.r!=moli.r) answer='n';
if(ensi.a!=moli.a) answer='n';
return answer;
}
/* Move to new file
*/
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
/* Move to new file
*/
void set_molecule_molindexes(molecule *m, int mi){
int i;
m[0].mi=mi;
for(i=0;i<m[0].nr;i++)
  {
printf("setting molecule %d residue indexed at %d\n",mi,i);
  set_residue_molindexes(&m[0].r[i],mi,i);
  }
}
/* Move to new file
*/
void set_assembly_molindexes(assembly *A){
int i;
for(i=0;i<A[0].nm;i++)
  {
  set_molecule_molindexes(&A[0].m[i][0],i);
  }
}
/* Move to new file
*/
void set_ensemble_molindexes(ensemble *E){
int i;
for(i=0;i<E[0].nm;i++)
  {
  set_molecule_molindexes(&E[0].m[i],i);
  }
}



/*
    This function sets molbonds at the residue level (r.rb).
*/
void set_molecule_residue_molbonds(molecule *m)
{
int ri,ai,bi,this_b;

for(ri=0;ri<m[0].nr;ri++)
    {
    /*
    Count the number of interresidue bonds.
    */
    m[0].r[ri].nrb=0;
    for(ai=0;ai<m[0].r[ri].na;ai++)
        {
        for(bi=0;bi<m[0].r[ri].a[ai].nmb;bi++)
            {
            if(m[0].r[ri].a[ai].mb[bi].s.r!=m[0].r[ri].a[ai].mb[bi].t.r) m[0].r[ri].nrb++;
            }
        }
    m[0].r[ri].rb=(molbond*)calloc(m[0].r[ri].nrb,sizeof(molbond));
    for(ai=0;ai<m[0].r[ri].na;ai++)
        {
        this_b=0;
        for(bi=0;bi<m[0].r[ri].a[ai].nmb;bi++)
            {
            if(m[0].r[ri].a[ai].mb[bi].s.r!=m[0].r[ri].a[ai].mb[bi].t.r)
                {
                m[0].r[ri].rb[this_b]=m[0].r[ri].a[ai].mb[bi];
                this_b++;
                }
            }
        if(this_b>m[0].r[ri].nrb){mywhine("this_b>m[0].r[ri].nrb in set_molecule_residue_molbonds.");}
        }
    }

return;
}

/*
    The  next two functions set the connection tree (m.rT) at the
        residue level (sets bonding between residues).  It uses r.rb
        so see also set_molecule_residue_molbonds().
*/
void set_molecule_residue_nodes_from_bonds(molecule *m)
{
int ri,finished;
char ensisame='y';
m[0].rT=(residue_node*)calloc(m[0].nr,sizeof(residue_node));
for(ri=0;ri<m[0].nr;ri++) 
    { 
    m[0].r[ri].mTi=-ri-1;
    m[0].rT[ri].isorigin='N';
    m[0].rT[ri].ID=copy_moli_to_ensi(m[0].r[ri].moli);
    if(ri>0)
        {
        ensisame=is_consistent_ensi_ensi(m[0].rT[ri].ID,m[0].rT[ri-1].ID);
/*printf("m-r-a ri=%d is %d - %d - %d \n",ri,m[0].rT[ri].ID.m,m[0].rT[ri].ID.r,m[0].rT[ri].ID.a);
printf("m-r-a ri-1=%d is %d - %d - %d \n",ri-1,m[0].rT[ri-1].ID.m,m[0].rT[ri-1].ID.r,m[0].rT[ri-1].ID.a);*/
        if(ensisame!='n')
            {
            mywhine("The residue molindexes are not unique in set_molecule_residue_nodes_from_bonds.\n");
            }
        }
    }
if(m[0].nr==1) 
    {
    m[0].r[0].mTi=0;
    return;
    }
m[0].rT[0].isorigin='Y';
finished = follow_molecule_residue_nodes_from_bonds(m, 0, &m[0].r[0]);
return;
}
int follow_molecule_residue_nodes_from_bonds(molecule *m, int iTree, residue *r)
{
char is_incoming='n';
int bi,ni,oi,finished=0;
residue *rt;
/* This function doesn't choose where to start.  It just follows. */

/* And, the value of mTi should be negative */
if(r[0].mTi>=0){mywhine("unexpected init of r[0].mTi in follow_molecule_residue_nodes_from_bonds.");}

/* 
0.  If this is not the end of a list, declare space for outgoing bonds 
*/
if(m[0].rT[iTree].nmbi<0){mywhine("Unexpected init of rT[].nmbi in follow_molecule_residue_nodes_from_bonds.");}
m[0].rT[iTree].nmbo=r[0].nrb-m[0].rT[iTree].nmbi;
/* 
    If this is the end of a list, go find any unseen or bail.
*/
if(m[0].rT[iTree].nmbo==0)
    {
    r[0].mTi=-(r[0].mTi+1);
    for(ni=0;ni<m[0].nr;ni++)
        {
        if(m[0].r[ni].mTi<0) 
           { 
           finished = follow_molecule_residue_nodes_from_bonds(m,ni,&m[0].r[ni]);
           }
        }
    return -1;
    }
if(finished==-1) { return -1; }
m[0].rT[iTree].mbo=(molbond**)calloc(m[0].rT[iTree].nmbo,sizeof(molbond*));
for(oi=0;oi<m[0].rT[iTree].nmbo;oi++)
    {
    m[0].rT[iTree].mbo[oi]=(molbond*)calloc(1,sizeof(molbond));
    }

/*
1.  Set each bond that is not an incoming bond as an outgoing bond 
*/
oi=0;
for(bi=0;bi<r[0].nrb;bi++)
    {
    /*
    If the two atoms are not part of the same molecule, complain, and don't bother.
    */
    if(r[0].rb[bi].s.m!=r[0].rb[bi].t.m) 
        {
        printf("Found bond between two molecules in follow_molecule_residue_nodes_from_bonds.\n");
        printf("Ignoring that bond.\n");
        m[0].rT[iTree].nmbo--;
        continue;
        }
    /*
    Check to see if the current bond is one of the incoming bonds.
    */
    is_incoming='n';
    for(ni=0;ni<m[0].rT[iTree].nmbi;ni++)
        {
        if(is_consistent_molbond_molbond(m[0].rT[iTree].mbi[ni][0],r[0].rb[bi])!='n') is_incoming='y';
        }
    /*
    If not, then set this as outgoing.
    */
    if(is_incoming=='n') 
        {
/*printf("m[0].rT[iTree].nmbo = %d\n",m[0].rT[iTree].nmbo);
*/
        m[0].rT[iTree].mbo[oi][0]=r[0].rb[bi];
        oi++;
        }
    if(oi>m[0].rT[iTree].nmbo){mywhine("overstepped outgoing bonds in follow_molecule_residue_nodes_from_bonds.");}
    }

/* 
2.  Set the current atom's state as seen.  Assumes a.rTi is initialized negative. 
*/
r[0].mTi = -(r[0].mTi + 1); /* make positive and add one*/

/*
3.  Identify the current atom as incoming to each outgoing atom. 
*/
for(oi=0;oi<m[0].rT[iTree].nmbo;oi++)
    {
    /* 
    get target residue
    */
    rt=&m[0].r[m[0].rT[iTree].mbo[oi][0].t.r]; 
/*printf("oi is %d ; target r=%d \n",oi, m[0].rT[iTree].mbo[oi][0].t.r);
*/
    /* 
    record location in contree 
    */
    if(rt[0].mTi>=0) 
        {
        printf("\nTarget residue mTi is nonnegative in follow_molecule_residue_nodes_from_bonds.\n");
        printf("This might be a problem.\n");
        bi=rt[0].mTi;
        } 
    else{bi=-(rt[0].mTi+1);}
    /* 
    identify current residue as incoming to the target atom
    */
    m[0].rT[bi].nmbi++; 
/*printf("rt[0].mTi is %d ; bi is %d ; m[0].rT[bi].nmbi is %d\n",rt[0].mTi,bi,m[0].rT[bi].nmbi);
*/
    if(m[0].rT[bi].nmbi==1) { m[0].rT[bi].mbi=(molbond**)calloc(1,sizeof(molbond*)); }
    else { m[0].rT[bi].mbi=(molbond**)realloc(m[0].rT[bi].mbi,m[0].rT[bi].nmbi*sizeof(molbond*));}
    m[0].rT[bi].mbi[m[0].rT[bi].nmbi-1]=(molbond*)calloc(1,sizeof(molbond));
    m[0].rT[bi].mbi[m[0].rT[bi].nmbi-1][0]=m[0].rT[iTree].mbo[oi][0];
    }

/* 
4.  Call each outgoing atom, in order, with this function 
*/
for(oi=0;oi<m[0].rT[iTree].nmbo;oi++)
    { 
    /* 
    get target residue
    */
    rt=&m[0].r[m[0].rT[iTree].mbo[oi][0].t.r]; 
    /* 
    record location in contree 
    */
    if(rt[0].mTi>=0) { bi=rt[0].mTi; } 
    else{bi=-(rt[0].mTi+1);}
    finished = follow_molecule_residue_nodes_from_bonds(m, bi, rt);
    }

return finished;
}


void set_residue_atom_nodes_from_bonds(residue *r)
{
int ai,finished=0;
char ensisame='y';
r[0].aT=(atom_node*)calloc(r[0].na,sizeof(atom_node));
for(ai=0;ai<r[0].na;ai++) 
    { 
    r[0].a[ai].rTi=-ai-1;
    r[0].aT[ai].isorigin='N';
    r[0].aT[ai].ID=copy_moli_to_ensi(r[0].a[ai].moli);
    if(ai>0)
        {
        ensisame=is_consistent_ensi_ensi(r[0].aT[ai].ID,r[0].aT[ai-1].ID);
        if(ensisame!='n')
            {
            mywhine("The atom molindexes are not unique in set_residue_atom_nodes_from_bonds.\n");
            }
        }
    }
if(r[0].na==1) 
    {
    r[0].a[0].rTi=0;
    return;
    }
r[0].aT[0].isorigin='Y';
finished = follow_residue_atom_nodes_from_bonds(r, 0, &r[0].a[0]);
return;
}
int follow_residue_atom_nodes_from_bonds(residue *r, int iTree, atom *a)
{
char is_incoming='n';
int bi,ni,oi,finished=0;
atom *at;
/* 
    This function doesn't choose where to start.  It just follows.
 
    At this point, the value of rTi should be negative or something is wrong.
*/
if(a[0].rTi>=0){mywhine("Unexpected init of a[0].rTi in follow_residue_atom_nodes_from_bonds.");}

/* 
0.  Declare space for outgoing bonds 
*/
if(r[0].aT[iTree].ni<0){mywhine("Unexpected init of aT[].ni in follow_residue_atom_nodes_from_bonds.");}
r[0].aT[iTree].no=a[0].nmb-r[0].aT[iTree].ni;
/*printf("r.N is %s ; a.N is %s ; a[0].rTi is %d ; no is %d ; ni is %d ; nmb is %d ",\
r[0].N,a[0].N,a[0].rTi,r[0].aT[iTree].no,r[0].aT[iTree].ni,a[0].nmb); */
/* 
    If this is the end of a list, go find any unseen or bail.
*/
if(r[0].aT[iTree].no<=0)
    {
/*printf("r.N is %s ; a.N is %s ; a[0].rTi is %d (changing to ",r[0].N,a[0].N,a[0].rTi); */
    a[0].rTi=-(a[0].rTi+1);
/*printf("%d)\n\tScanning rTi's:",a[0].rTi); */
    for(ni=0;ni<r[0].na;ni++)
        {
/*printf("\t\t a[%d].N is %s ; a[0].rTi is %d \n",ni,r[0].a[ni].N,r[0].a[ni].rTi); */
        if(r[0].a[ni].rTi<0) 
           { 
/*printf("About to follow residue for r[0].a[ni].rTi=%d",r[0].a[ni].rTi); */
           finished = follow_residue_atom_nodes_from_bonds(r,ni,&r[0].a[ni]);
           }
        if(finished==-1) break;
        }
    if (finished == -1) { return -1;}
    if (ni==r[0].na) { return -1;}
    else { return -2;}
    } 
/* 
    At this point, the value of rTi should be negative or something is wrong.
*/
if(a[0].rTi>=0){mywhine("Unexpected init of a[0].rTi in follow_residue_atom_nodes_from_bonds.");}

r[0].aT[iTree].o=(ensindex**)calloc(r[0].aT[iTree].no,sizeof(ensindex*));
for(oi=0;oi<r[0].aT[iTree].no;oi++)
    {
    r[0].aT[iTree].o[oi]=(ensindex*)calloc(1,sizeof(ensindex));
    }

/* 
1.  Set each bond that is not an incoming bond as an outgoing bond 
*/
oi=0;
for(bi=0;bi<a[0].nmb;bi++)
    {
    /*
    If the two atoms are not part of the same residue, don't bother.
    */
/*printf("a.N is %s ; a[0].mb[bi].s.r=%d a[0].mb[bi].t.r=%d\n",a[0].N,a[0].mb[bi].s.r,a[0].mb[bi].t.r); */
    if(a[0].mb[bi].s.r!=a[0].mb[bi].t.r)
        {
        r[0].aT[iTree].no--;
/*printf("....skipping.   r[0].aT[iTree].no=%d \n",r[0].aT[iTree].no); */
        continue;
        }
    /*
    Check to see if the current bond is one of the incoming bonds.
    */
    is_incoming='n';
    for(ni=0;ni<r[0].aT[iTree].ni;ni++)
        {
/*printf("This is the incoming check\n"); */
/*printf("\tThe result of is_consistent_ensi_moli is %c\n",\
  is_consistent_ensi_moli(r[0].aT[iTree].i[ni][0],a[0].mb[bi].t)); */
/*printf("\tThe values are:\n"); */
/*printf("\t\tr[0].aT[iTree].i[ni][0] -- i,m,r,a = %d %d %d %d \n",\
  r[0].aT[iTree].i[ni][0].i,r[0].aT[iTree].i[ni][0].m,r[0].aT[iTree].i[ni][0].r,r[0].aT[iTree].i[ni][0].a); */
/*printf("\t\ta[0].mb[bi].t -- i,m,r,a = %d %d %d %d \n",\
  a[0].mb[bi].t.i,a[0].mb[bi].t.m,a[0].mb[bi].t.r,a[0].mb[bi].t.a); */
        if(is_consistent_ensi_moli(r[0].aT[iTree].i[ni][0],a[0].mb[bi].t)!='n') is_incoming='y';
        }
/*printf("past incoming check, is_incoming=%c  oi is %d\n",is_incoming,oi); */
    /*
    If not, then set this as outgoing.
    */
    if(is_incoming=='n') 
        {
/*printf("The target ID is a.N=%s\n",r[0].a[a[0].mb[bi].t.a].N); */
        r[0].aT[iTree].o[oi][0]=copy_moli_to_ensi(a[0].mb[bi].t);
        oi++;
        }
    if(oi>r[0].aT[iTree].no){mywhine("overstepped outgoing bonds in follow_residue_atom_nodes_from_bonds.");}
    }

/* 
2.  Set the current atom's state as seen.  Assumes a.rTi is initialized negative. 
*/
a[0].rTi = -(a[0].rTi + 1); /* make positive and add one*/

/*
3.  Identify the current atom as incoming to each outgoing atom. 
*/
for(oi=0;oi<r[0].aT[iTree].no;oi++)
    {
    /* 
    get target atom 
    */
    at=&r[0].a[r[0].aT[iTree].o[oi][0].a]; 
    /* 
    record location in contree 
    */
    if(at[0].rTi>=0) { continue; } 
    else{bi=-(at[0].rTi+1);}
    /* 
    identify current atom as incoming to the target atom
    */
    r[0].aT[bi].ni++; 
    if(r[0].aT[bi].ni==1) { r[0].aT[bi].i=(ensindex**)calloc(1,sizeof(ensindex*)); }
    else { r[0].aT[bi].i=(ensindex**)realloc(r[0].aT[bi].i,r[0].aT[bi].ni*sizeof(ensindex*));}
    r[0].aT[bi].i[r[0].aT[bi].ni-1]=(ensindex*)calloc(1,sizeof(ensindex));
    r[0].aT[bi].i[r[0].aT[bi].ni-1][0]=r[0].aT[iTree].ID;
    /*
    */
    }

/* 
4.  Call each outgoing atom, in order, with this function 
*/
for(oi=0;oi<r[0].aT[iTree].no;oi++)
    { 
    /* 
    get target atom 
    */
    at=&r[0].a[r[0].aT[iTree].o[oi][0].a]; 
    /* 
    record location in contree 
    */
    if(at[0].rTi>=0) { continue; } 
    else{bi=-(at[0].rTi+1);}
/*printf("About to call follow_residue_atom_nodes_from_bonds for:\n"); */
/*printf("\tr.N=%s ; at.N=%s ; bi = %d ; at[0].rTi = %d\n",r[0].N,at[0].N,bi,at[0].rTi); */
    finished = follow_residue_atom_nodes_from_bonds(r, bi, at);
    }

return finished;
}

int  follow_molecule_atom_nodes_from_bonds(molecule *m, int iTree, atom *a)
{
char is_incoming='n';
int bi,ni,oi;
residue *rt;
atom *at;
/* This function doesn't choose where to start.  It just follows. */

/* And, the value of rTi should be negative */
if(a[0].mTi>=0){mywhine("unexpected init of a[0].rTi in follow_molecule_atom_nodes_from_bonds.");}

/* 
0.  Declare space for outgoing bonds 
*/
if(m[0].aT[iTree].ni<0){mywhine("Unexpected init of aT[].ni in follow_molecule_atom_nodes_from_bonds");}
m[0].aT[iTree].no=a[0].nmb-m[0].aT[iTree].ni;
m[0].aT[iTree].o=(ensindex**)calloc(m[0].aT[iTree].no,sizeof(ensindex*));
for(oi=0;oi<m[0].aT[iTree].no;oi++)
    {
    m[0].aT[iTree].o[oi]=(ensindex*)calloc(1,sizeof(ensindex));
    }

/* 
1.  Set each bond that is not an incoming bond as an outgoing bond 
*/
oi=0;
for(bi=0;bi<a[0].nmb;bi++)
    {
    /*
    If the two atoms are not part of the same molecule, complain, and don't bother.
    */
    if(a[0].mb[bi].s.m!=a[0].mb[bi].t.m) 
        {
        printf("Found bond between two molecules in follow_molecule_atom_nodes_from_bonds.\n");
        printf("Ignoring that bond.\n");
        m[0].aT[iTree].no--;
        continue;
        }
    /*
    Check to see if the current bond is one of the incoming bonds.
    */
    is_incoming='n';
    for(ni=0;ni<m[0].aT[iTree].ni;ni++)
        {
        if(is_consistent_ensi_moli(m[0].aT[iTree].i[ni][0],a[0].mb[bi].t)!='n') is_incoming='y';
        }
    /*
    If not, then set this as outgoing.
    */
    if(is_incoming=='n') 
        {
        m[0].aT[iTree].o[oi][0]=copy_moli_to_ensi(a[0].mb[bi].t);
        oi++;
        }
    if(oi>m[0].aT[iTree].no){mywhine("overstepped outgoing bonds in follow_molecule_atom_nodes_from_bonds.");}
    }

/* 
2.  Set the current atom's state as seen.  Assumes a.rTi is initialized negative. 
*/
a[0].rTi = -(a[0].mTi + 1); /* make positive and add one*/

/*
3.  Identify the current atom as incoming to each outgoing atom. 
*/
for(oi=0;oi<m[0].aT[iTree].no;oi++)
    {
    /* 
    get target atom 
    */
    at=&m[0].r[m[0].aT[iTree].o[oi][0].r].a[m[0].aT[iTree].o[oi][0].a]; 
    /* 
    record location in contree 
    */
    if(at[0].mTi>=0) 
        {
        printf("\nTarget atom mTi is nonnegative in follow_molecule_atom_nodes_from_bonds.\n");
        printf("This might be a problem.\n");
        bi=at[0].mTi;
        } 
    else{bi=-(at[0].mTi+1);}
    /* 
    identify current atom as incoming to the target atom
    */
    m[0].aT[bi].ni++; 
    if(m[0].aT[bi].ni==1) { m[0].aT[bi].i=(ensindex**)calloc(1,sizeof(ensindex*)); }
    else { m[0].aT[bi].i=(ensindex**)realloc(m[0].aT[bi].i,m[0].aT[bi].ni*sizeof(ensindex*));}
    m[0].aT[bi].i[m[0].aT[bi].ni-1]=(ensindex*)calloc(1,sizeof(ensindex));
    m[0].aT[bi].i[m[0].aT[bi].ni-1][0]=m[0].aT[iTree].ID;
    /*
    */
    }

/* 
4.  Call each outgoing atom, in order, with this function 
*/
for(oi=0;oi<m[0].aT[iTree].no;oi++)
    { 
    /* 
    get target atom 
    */
    rt=&m[0].r[m[0].aT[iTree].o[oi][0].r]; 
    at=&rt[0].a[m[0].aT[iTree].o[oi][0].a]; 
    /* 
    record location in contree 
    */
    if(at[0].mTi>=0) { bi=at[0].mTi; } 
    else{bi=-(at[0].mTi+1);}
    follow_molecule_atom_nodes_from_bonds(m, bi, at);
    }

return 0;
} 
