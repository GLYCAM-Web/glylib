/** \file bonding_utilities.c
*   Purpose:  Functions associated with bonding. 
*   Begin 20110417 BLFoley
*
*   In general, the idea is to set bonding within residues, then translate
*   that up to the molecule level.
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

char is_consistent_ensi_moli(ensindex ensi, molindex moli){
char answer='y';
if(ensi.i!=moli.i) answer='i';
if(ensi.m!=moli.m) answer='n';
if(ensi.r!=moli.r) answer='n';
if(ensi.a!=moli.a) answer='n';
return answer;
}

int follow_residue_atom_bonds_for_contree(residue *r);
int follow_residue_atom_bonds_for_contree(residue *r, int iTree, atom *a){
char is_incoming='n';
int bi,ni,oi;
/* This function doesn't choose where to start.  It just follows. */

/* 0.  Declare space for outgouing bonds */
if(r[0].aT[iTree].ni<0){mywhine("Unexpected init of aT[].ni in residue_atom_bonds.");}
r[0].aT[iTree].no=a[0].nmb-r[0].aT[iTree].ni;
r[0].aT[iTree].o=(ensindex**)calloc(r[0].aT[iTree].no,sizeof(ensindex*));
for(oi=0;oi<r[0].aT[iTree].no;oi++)
    {
    r[0].aT[iTree].o[oi]=(ensindex*)calloc(1,sizeof(ensindex));
    }

/* 1.  Set each bond that is not an incoming bond as an outgoing bond */
oi=0;
for(bi=0;bi<a[0].nmb;bi++)
    {
    if(a[0].mb[bi].s.r!=a[0].mb[bi].t.r) continue;
    is_incoming='n';
    for(ni=0;ni<r[0].aT[iTree].ni;ni++)
        {
        if(is_consistent_ensi_moli(r[0].aT[iTree].i[ni][0],a[0].mb[bi])!='n') is_incoming='y';
        }
    if(is_incoming=='n') 
        {
        r[0].aT[iTree].
        }
    }

/* 2.  Set the current atom's state as seen. */
/* 3.  Call each outgoing atom, in order, with this function */

return NextNumber;
}
