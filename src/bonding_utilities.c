/** \file bonding_utilities.c
*   Purpose:  Functions associated with bonding. 
*   Begin 20110417 BLFoley
*
*/
void set_molecule_atom_nodes_from_bonds(molecule *m, int iTree, atom *a);
void set_residue_atom_nodes_from_bonds(residue *r, int iTree, atom *a);
void set_molecule_residue_molbonds(molecule *m);
void set_molecule_residue_nodes_from_bonds(molecule *m, int iTree, residue *r);
ensindex copy_moli_to_ensi(molindex moli);
char is_consistent_ensi_moli(ensindex ensi, molindex moli);

ensindex copy_moli_to_ensi(molindex moli){
ensindex ensi;
ensi.E=ensi.A=-1;
ensi.i=moli.i;
ensi.m=moli.m;
ensi.r=moli.r;
ensi.a=moli.a;
return ensi;
}

char is_consistent_ensi_moli(ensindex ensi, molindex moli){
char answer='y';
if(ensi.i!=moli.i) answer='i';
if(ensi.m!=moli.m) answer='n';
if(ensi.r!=moli.r) answer='n';
if(ensi.a!=moli.a) answer='n';
return answer;
}

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

void set_molecule_residue_nodes_from_bonds(molecule *m, int iTree, residue *r)
{
char is_incoming='n';
int bi,ni,oi;
atom *at;
residue *rt;
/* This function doesn't choose where to start.  It just follows. */

/* And, the value of mTi should be negative */
if(r[0].mTi>=0){mywhine("unexpected init of r[0].mTi in set_molecule_residue_nodes_from_bonds.");}

/* 
0.  Declare space for outgoing bonds 
*/
if(m[0].rT[iTree].ni<0){mywhine("Unexpected init of rT[].ni in set_molecule_residue_nodes_from_bonds.");}
m[0].rT[iTree].no=r[0].nrb-m[0].rT[iTree].ni;
m[0].rT[iTree].o=(ensindex**)calloc(m[0].rT[iTree].no,sizeof(ensindex*));
for(oi=0;oi<m[0].rT[iTree].no;oi++)
    {
    m[0].rT[iTree].o[oi]=(ensindex*)calloc(1,sizeof(ensindex));
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
        printf("Found bond between two molecules in set_molecule_residue_nodes_from_bonds.\n");
        printf("Ignorning that bond.\n");
        continue;
        }
    /*
    Check to see if the current bond is one of the incoming bonds.
    */
    is_incoming='n';
    for(ni=0;ni<m[0].rT[iTree].ni;ni++)
        {
        if(is_consistent_ensi_moli(m[0].rT[iTree].i[ni][0],r[0].rb[bi])!='n') is_incoming='y';
        }
    /*
    If not, then set this as outgoing.
    */
    if(is_incoming=='n') 
        {
        m[0].rT[iTree].o[oi][0]=copy_moli_to_ensi(r[0].rb[bi]);
        oi++;
        }
    if(oi>m[0].rT[iTree].no){mywhine("overstepped outgoing bonds in set_molecule_residue_nodes_from_bonds.");}
    }

/* 
2.  Set the current atom's state as seen.  Assumes a.rTi is initialized negative. 
*/
r[0].mTi = -r[0].mTi + 1; /* make positive and add one*/

/*
3.  Identify the current atom as incoming to each outgoing atom. 
*/
for(oi=0;oi<r[0].aT[iTree].no;oi++)
    {
    /*
    If the two atoms are not part of the same molecule, don't bother.
    */
    if(r[0].rb[bi].s.m!=r[0].rb[bi].t.m) continue; 
    /* 
    get target atom and residue
    */
    rt=&m[0].r[r[0].mT[iTree].o[oi][0].r]; 
    at=&rt[0].a[r[0].mT[iTree].o[oi][0].a]; 
    /* 
    record location in contree 
    */
    if(rt[0].mTi>=0) 
        {
        printf("\nTarget residue mTi is nonnegative in set_molecule_residue_nodes_from_bonds.\n");
        printf("This might be a problem.\n");
        bi=rt[0].mTi;
        } 
    else{bi=-rt[0].mTi+1;}
    /* 
    identify current residue as incoming to the target atom
    */
    m[0].rT[bi].ni++; 
    if(m[0].rT[bi].ni==1) { m[0].rT[bi].i=(ensindex**)calloc(1,sizeof(ensindex*)); }
    else { m[0].rT[bi].i=(ensindex**)realloc(m[0].rT[bi].i,m[0].rT[bi].ni*sizeof(ensindex*));}
    m[0].rT[bi].i[m[0].rT[bi].ni-1]=(ensindex*)calloc(1,sizeof(ensindex));
    m[0].rT[bi].i[m[0].rT[bi].ni-1][0]=m[0].rT[iTree].ID;
    }

/* 

START HERE

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
    if(at[0].rTi>=0) { bi=at[0].rTi; } 
    else{bi=-at[0].rTi+1;}
    follow_molecule_residue_nodes_from_bonds(r, bi, at);
    }

return;
}

void set_residue_atom_nodes_from_bonds(residue *r, int iTree, atom *a)
{
char is_incoming='n';
int bi,ni,oi;
atom *at;
/* This function doesn't choose where to start.  It just follows. */

/* And, the value of rTi should be negative */
if(a[0].rTi>=0){mywhine("unexpected init of a[0].rTi in follow_residue_atom_bonds_for_nodes..");}

/* 
0.  Declare space for outgoing bonds 
*/
if(r[0].aT[iTree].ni<0){mywhine("Unexpected init of aT[].ni in follow_residue_atom_bonds_for_nodes..");}
r[0].aT[iTree].no=a[0].nmb-r[0].aT[iTree].ni;
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
    if(a[0].mb[bi].s.r!=a[0].mb[bi].t.r) continue;
    /*
    Check to see if the current bond is one of the incoming bonds.
    */
    is_incoming='n';
    for(ni=0;ni<r[0].aT[iTree].ni;ni++)
        {
        if(is_consistent_ensi_moli(r[0].aT[iTree].i[ni][0],a[0].mb[bi])!='n') is_incoming='y';
        }
    /*
    If not, then set this as outgoing.
    */
    if(is_incoming=='n') 
        {
        r[0].aT[iTree].o[oi][0]=copy_moli_to_ensi(a[0].mb[bi]);
        oi++;
        }
    if(oi>r[0].aT[iTree].no){mywhine("overstepped outgoing bonds in follow_residue_atom_bonds_for_nodes.");}
    }

/* 
2.  Set the current atom's state as seen.  Assumes a.rTi is initialized negative. 
*/
a[0].rTi = -a[0].rTi + 1; /* make positive and add one*/

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
    if(at[0].rTi>=0) 
        {
        printf("\nTarget atom rTi is nonnegative in follow_residue_atom_bonds_for_nodes.\n");
        printf("This might be a problem.\n");
        bi=at[0].rTi;
        } 
    else{bi=-at[0].rTi+1;}
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
    if(at[0].rTi>=0) { bi=at[0].rTi; } 
    else{bi=-at[0].rTi+1;}
    follow_residue_atom_bonds_for_nodes(r, bi, at);
    }

return;
}

void set_molecule_atom_nodes_from_bonds(molecule *m, int iTree, atom *a)
{
char is_incoming='n';
int bi,ni,oi;
atom *at;
/* This function doesn't choose where to start.  It just follows. */

/* And, the value of rTi should be negative */
if(a[0].mTi>=0){mywhine("unexpected init of a[0].rTi in set_molecule_atom_nodes_from_bonds.");}

/* 
0.  Declare space for outgoing bonds 
*/
if(m[0].aT[iTree].ni<0){mywhine("Unexpected init of aT[].ni in set_molecule_atom_nodes_from_bonds");}
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
        printf("Found bond between two molecules in set_molecule_atom_nodes_from_bonds.\n");
        printf("Ignorning that bond.\n");
        continue;
        }
    /*
    Check to see if the current bond is one of the incoming bonds.
    */
    is_incoming='n';
    for(ni=0;ni<m[0].aT[iTree].ni;ni++)
        {
        if(is_consistent_ensi_moli(m[0].aT[iTree].i[ni][0],a[0].mb[bi])!='n') is_incoming='y';
        }
    /*
    If not, then set this as outgoing.
    */
    if(is_incoming=='n') 
        {
        m[0].aT[iTree].o[oi][0]=copy_moli_to_ensi(a[0].mb[bi]);
        oi++;
        }
    if(oi>m[0].aT[iTree].no){mywhine("overstepped outgoing bonds in follow_residue_atom_bonds_for_nodes.");}
    }

/* 
2.  Set the current atom's state as seen.  Assumes a.rTi is initialized negative. 
*/
a[0].rTi = -a[0].mTi + 1; /* make positive and add one*/

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
        printf("\nTarget atom mTi is nonnegative in follow_molecule_atom_bonds_for_nodes.\n");
        printf("This might be a problem.\n");
        bi=at[0].mTi;
        } 
    else{bi=-at[0].mTi+1;}
    /* 
    identify current atom as incoming to the target atom
    */
    m[0].aT[bi].ni++; 
    if(m[0].aT[bi].ni==1) { m[0].aT[bi].i=(ensindex**)calloc(1,sizeof(ensindex*)); }
    else { m[0].aT[bi].i=(ensindex**)realloc(m[0].aT[bi].i,m[0].aT[bi].ni*sizeof(ensindex*));}
    m[0].aT[bi].i[m[0].aT[bi].ni-1]=(ensindex*)calloc(1,sizeof(ensindex));
    m[0].aT[bi].i[m[0].aT[bi].ni-1][0]=copy_moli_to_ensi(m[0].aT[iTree].ID);
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
    at=&m[0].a[m[0].aT[iTree].o[oi][0].a]; 
    /* 
    record location in contree 
    */
    if(at[0].mTi>=0) { bi=at[0].mTi; } 
    else{bi=-at[0].mTi+1;}
    follow_residue_atom_bonds_for_nodes(m, bi, at);
    }

return;
}



