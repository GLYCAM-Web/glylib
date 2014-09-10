/** \file bonding_utilities.c
* \addtogroup BONDING
*   Purpose:  Functions associated with bonding. 
*   Begin 20110417 BLFoley
*
*   Return values for follow functions:
*        0 : Made it all the way to the end.  All ok.
*       -1 : Made it to an appropriate end in the middle.  Don't continue.
*       -2 : Something else happened, but probably is ok.
* 
*/
#include <mylib.h>
#include <molecules.h>

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
    this_b=0;
    for(ai=0;ai<m[0].r[ri].na;ai++)
        {
/*printf("This atom is %s\n",m[0].r[ri].a[ai].N);*/
        for(bi=0;bi<m[0].r[ri].a[ai].nmb;bi++)
            {
/*printf("\tm[0].r[%d].a[%d].mb[%d].s.r = %d\n",ri,ai,bi,m[0].r[ri].a[ai].mb[bi].s.r);*/
/*printf("\tm[0].r[%d].a[%d].mb[%d].t.r = %d\n",ri,ai,bi,m[0].r[ri].a[ai].mb[bi].t.r);*/
            if(m[0].r[ri].a[ai].mb[bi].s.r!=m[0].r[ri].a[ai].mb[bi].t.r)
                {
                m[0].r[ri].rb[this_b]=m[0].r[ri].a[ai].mb[bi];
/*printf("\t this_b is %d\n",this_b);*/
/*printf("\tfound molbond for %s from %d-%d-%d-%d to %d-%d-%d-%d\n",m[0].r[ri].N,\
m[0].r[ri].a[ai].mb[bi].s.i,m[0].r[ri].a[ai].mb[bi].s.m,m[0].r[ri].a[ai].mb[bi].s.r,m[0].r[ri].a[ai].mb[bi].s.a,\
m[0].r[ri].a[ai].mb[bi].t.i,m[0].r[ri].a[ai].mb[bi].t.m,m[0].r[ri].a[ai].mb[bi].t.r,m[0].r[ri].a[ai].mb[bi].t.a);*/
/*printf("\tsetting it to residue level as %d-%d-%d-%d to %d-%d-%d-%d\n",\
m[0].r[ri].rb[this_b].s.i,m[0].r[ri].rb[this_b].s.m,m[0].r[ri].rb[this_b].s.r,m[0].r[ri].rb[this_b].s.a,\
m[0].r[ri].rb[this_b].t.i,m[0].r[ri].rb[this_b].t.m,m[0].r[ri].rb[this_b].t.r,m[0].r[ri].rb[this_b].t.a);*/
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
int ri;
char ensisame='y';
m[0].rT=(residue_node*)calloc(m[0].nr,sizeof(residue_node));
for(ri=0;ri<m[0].nr;ri++) 
    { 
    m[0].r[ri].mTi=-ri-1;
    m[0].rT[ri].isorigin='N';
    m[0].rT[ri].ID=copy_moli_to_ensi(m[0].r[ri].moli);
    m[0].rT[ri].ID.i=ri;
    if(ri>0)
        {
        ensisame=is_consistent_ensi_ensi(m[0].rT[ri].ID,m[0].rT[ri-1].ID);
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
follow_molecule_residue_nodes_from_bonds(m, 0, &m[0].r[0]);
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
if(m[0].rT[iTree].nmbi<0)
  {
  mywhine("Unexpected init of rT[].nmbi in follow_molecule_residue_nodes_from_bonds.");
  }
m[0].rT[iTree].nmbo=r[0].nrb-m[0].rT[iTree].nmbi;
/* 
    If this is the end of a list, go find any unseen or bail.
*/
if(m[0].rT[iTree].nmbo<=0)
    {
    r[0].mTi=-(r[0].mTi+1);
    for(ni=0;ni<m[0].nr;ni++)
        {
        if(m[0].r[ni].mTi<0) 
           { 
           finished = follow_molecule_residue_nodes_from_bonds(m,ni,&m[0].r[ni]);
           }
        if(finished==-1) break;
        }
    if (finished == -1) { return -1;}
    if (ni==m[0].nr) { return -1;}
    else { return -2;}
    }

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
        if(is_consistent_molbond_molbond_inverse(m[0].rT[iTree].mbi[ni][0],r[0].rb[bi])!='n') is_incoming='y';
        }
    /*
    If not, then set this as outgoing.
    */
    if(is_incoming=='n') 
        {
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
    m[0].rT[iTree].mbo[oi][0].i=bi;
    /* 
    identify current residue as incoming to the target atom
    */
    m[0].rT[bi].nmbi++; 
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
    if(rt[0].mTi>=0) { continue; } 
    else{bi=-(rt[0].mTi+1);}
    finished = follow_molecule_residue_nodes_from_bonds(m, bi, rt);
    }

return finished;
}


void set_residue_atom_nodes_from_bonds(residue *r)
{
int ai;
char ensisame='y';
r[0].aT=(atom_node*)calloc(r[0].na,sizeof(atom_node));
for(ai=0;ai<r[0].na;ai++) 
    { 
    r[0].a[ai].rTi=-ai-1;
    r[0].aT[ai].isorigin='N';
    r[0].aT[ai].ID=copy_moli_to_ensi(r[0].a[ai].moli);
    r[0].aT[ai].ID.i=ai;
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
follow_residue_atom_nodes_from_bonds(r, 0, &r[0].a[0]);
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
/* 
    If this is the end of a list, go find any unseen or bail.
*/
if(r[0].aT[iTree].no<=0)
    {
    a[0].rTi=-(a[0].rTi+1);
    for(ni=0;ni<r[0].na;ni++)
        {
        if(r[0].a[ni].rTi<0) 
           { 
           finished = follow_residue_atom_nodes_from_bonds(r,ni,&r[0].a[ni]);
           }
        if(finished==-1) break;
        }
    if (finished == -1) { return -1;}
    if (ni==r[0].na) { return -1;}
    else { return -2;}
    } 

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
    if(a[0].mb[bi].s.r!=a[0].mb[bi].t.r)
        {
        r[0].aT[iTree].no--;
        continue;
        }
    /*
    Check to see if the current bond is one of the incoming bonds.
    */
    is_incoming='n';
    for(ni=0;ni<r[0].aT[iTree].ni;ni++)
        {
        if(is_consistent_ensi_moli(r[0].aT[iTree].i[ni][0],a[0].mb[bi].t)!='n') is_incoming='y';
        }
    /*
    If not, then set this as outgoing.
    */
    if(is_incoming=='n') 
        {
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
    r[0].aT[iTree].o[oi][0].i=bi;
    /* 
    identify current atom as incoming to the target atom
    */
    r[0].aT[bi].ni++; 
    if(r[0].aT[bi].ni==1) { r[0].aT[bi].i=(ensindex**)calloc(1,sizeof(ensindex*)); }
    else { r[0].aT[bi].i=(ensindex**)realloc(r[0].aT[bi].i,r[0].aT[bi].ni*sizeof(ensindex*));}
    r[0].aT[bi].i[r[0].aT[bi].ni-1]=(ensindex*)calloc(1,sizeof(ensindex));
    r[0].aT[bi].i[r[0].aT[bi].ni-1][0]=r[0].aT[iTree].ID;
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
    finished = follow_residue_atom_nodes_from_bonds(r, bi, at);
    }

return finished;
}

void set_molecule_atom_nodes_from_bonds(molecule *m)
{
int ri,ai,na=0,aai;
char ensisame='y'; 

for(ri=0;ri<m[0].nr;ri++)
    {
    na+=m[0].r[ri].na;
    }
if((m[0].na>0)&&(m[0].na!=na))
    {
    printf("\nWarning:  m.na != sum of all r.na in set_molecule_atom_nodes_from_bonds\n");
    printf("            Setting those to be equal.  Hope that\'s ok.\n");
    }
m[0].na=na;
m[0].aT=(atom_node*)calloc(m[0].na,sizeof(atom_node));
aai=0;
for(ri=0;ri<m[0].nr;ri++) 
    { 
for(ai=0;ai<m[0].r[ri].na;ai++) 
    { 
    m[0].r[ri].a[ai].mTi=-aai-1;
    m[0].aT[aai].isorigin='N';
    m[0].aT[aai].ID=copy_moli_to_ensi(m[0].r[ri].a[ai].moli);
    m[0].aT[aai].ID.i=aai;
    /*  The following isn't a comprehensive check.  Just a spot-check. */
    if(aai>0)
        {
        ensisame=is_consistent_ensi_ensi(m[0].aT[aai].ID,m[0].aT[aai-1].ID);
        if(ensisame!='n')
            {
            mywhine("The atom molindexes are not unique in set_molecule_atom_nodes_from_bonds.\n");
            }
        }
    aai++;
    }
    }
if(m[0].na==1) 
    {
    m[0].r[0].a[0].mTi=0;
    return;
    }
m[0].aT[0].isorigin='Y';
follow_molecule_atom_nodes_from_bonds(m, 0, &m[0].r[0].a[0]);
return;
}

int  follow_molecule_atom_nodes_from_bonds(molecule *m, int iTree, atom *a)
{
char is_incoming='n';
int bi,ni,oi,finished=0;
atom *at;
/* 
    This function doesn't choose where to start.  It just follows. 

    At this point, the value of rTi should be negative 
*/
if(a[0].mTi>=0){mywhine("unexpected init of a[0].mTi in follow_molecule_atom_nodes_from_bonds.");}

/* 
0.  Declare space for outgoing bonds 
*/
if(m[0].aT[iTree].ni<0){mywhine("Unexpected init of aT[].ni in follow_molecule_atom_nodes_from_bonds");}
m[0].aT[iTree].no=a[0].nmb-m[0].aT[iTree].ni;
/* 
    If this is the end of a list, go find any unseen or bail.
*/

if(m[0].aT[iTree].no<=0)
    {
    a[0].mTi=-(a[0].mTi+1);
    for(ni=0;ni<m[0].na;ni++)
        {
        at=&m[0].r[m[0].aT[ni].ID.r].a[m[0].aT[ni].ID.a];
        if(at[0].mTi<0) 
           { 
           finished = follow_molecule_atom_nodes_from_bonds(m,ni,at);
           }
        if(finished==-1) break;
        }
    if (finished == -1) { return -1;}
    if (ni==m[0].na) { return -1;}
    else { return -2;}
    } 

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
a[0].mTi = -(a[0].mTi + 1); /* make positive and add one*/

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
    if(at[0].mTi>=0) { continue; } 
    else{bi=-(at[0].mTi+1);}
    m[0].aT[iTree].o[oi][0].i=bi;
    /* 
    identify current atom as incoming to the target atom
    */
    m[0].aT[bi].ni++; 
    if(m[0].aT[bi].ni==1) { m[0].aT[bi].i=(ensindex**)calloc(1,sizeof(ensindex*)); }
    else { m[0].aT[bi].i=(ensindex**)realloc(m[0].aT[bi].i,m[0].aT[bi].ni*sizeof(ensindex*));}
    m[0].aT[bi].i[m[0].aT[bi].ni-1]=(ensindex*)calloc(1,sizeof(ensindex));
    m[0].aT[bi].i[m[0].aT[bi].ni-1][0]=m[0].aT[iTree].ID;
    }

/* 
4.  Call each outgoing atom, in order, with this function 
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
    if(at[0].mTi>=0) { continue; } 
    else{bi=-(at[0].mTi+1);}
    finished = follow_molecule_atom_nodes_from_bonds(m, bi, at);
    }

return finished;
}
