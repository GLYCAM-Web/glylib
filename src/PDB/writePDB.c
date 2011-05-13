/* Author: Michael Tessier
	Extended and updated by BLFoley */
#include <mylib.h>
#include <PDB.h>
#include <gly_codeutils.h>
#define TRUE  1
#define FALSE 0

fileslurp get_assembly_PDB_ATOM_lines(assembly *A,char isource,int savei,char raltname,int xs)
{
int mi,ri,ai,ainit,rinit,ntot,Li,nModels,Mi,this_x,nAlt,iFM;
fileslurp FA,*FM; 
char *aptr,*rptr,*tmp;
/* allocate savei space if not -1 and if not already allocated */
if(savei!=-1)
    {
    for(mi=0;mi<A[0].nm;mi++)
        {
        for(ri=0;ri<A[0].m[mi][0].nr;ri++)
            {
            if(A[0].m[mi][0].r[ri].ni==0) 
                {
                A[0].m[mi][0].r[ri].ni=savei+1;
                A[0].m[mi][0].r[ri].i=(int*)calloc(savei+1,sizeof(int));
                }
            else if (A[0].m[mi][0].r[ri].ni<savei+1)
                {
                A[0].m[mi][0].r[ri].ni=savei+1;
                A[0].m[mi][0].r[ri].i=\
                   (int*)realloc(A[0].m[mi][0].r[ri].i,(savei+1)*sizeof(int));
                }
            for(ai=0;ai<A[0].m[mi][0].r[ri].na;ai++)
                {
                if(A[0].m[mi][0].r[ri].a[ai].ni==0) 
                    {
                    A[0].m[mi][0].r[ri].a[ai].ni=savei+1;
                    A[0].m[mi][0].r[ri].a[ai].i=(int*)calloc(savei+1,sizeof(int));
                    }
                else if (A[0].m[mi][0].r[ri].a[ai].ni<savei+1)
                    {
                    A[0].m[mi][0].r[ri].a[ai].ni=savei+1;
                    A[0].m[mi][0].r[ri].a[ai].i=\
                       (int*)realloc(A[0].m[mi][0].r[ri].a[ai].i,(savei+1)*sizeof(int));
                    }
        } } }
    }

/* 
    As a quick error check, see if the first residue's altname exists, if requested. 
*/
if((tolower(raltname)=='y')&&(A[0].m[0][0].r[0].altname==NULL))
    {
    printf("\nYou requested an assembly PDB print using alternate residue names.\n");
    printf("But the first one has not been set.  Not checking others.\n");
    printf("This might cause trouble. (in get_assembly_ATOM_lines)\n");
    }

/* 
  Figure out how many slurps we will need based on the coords to be printed.
*/
if(xs==-3)  /* print main plus all alternates */  
  { 
  nAlt=A[0].m[0][0].r[0].a[0].nalt; /* for error-checking, grab number alts for first atom */
  nModels=nAlt+1; 
  this_x=-1;
  }
else if (xs==-2) /* print only all alternates */
  {
  nAlt=A[0].m[0][0].r[0].a[0].nalt; /* for error-checking, grab number alts for first atom */
  nModels=nAlt; 
  this_x=0; 
  }
else if (xs==-1) /* print only the main set of coordinates */
  {
  nAlt=0; 
  nModels=1; 
  this_x=-1; 
  }
else /* print only one alternate set */
  {
  nAlt=xs+1; 
  nModels=1; 
  this_x=xs; 
  }

FM=(fileslurp*)calloc((nModels*A[0].nm),sizeof(fileslurp));
if(tolower(isource)=='n') ainit=rinit=-1;
else ainit=rinit=1;

ntot=0; /* the total number of lines returned by get_molecule_PDB_ATOM_lines */

for(Mi=0;Mi<nModels;Mi++)
    {
    iFM=Mi*A[0].nm;
for(mi=0;mi<A[0].nm;mi++)
    {
    /*
        Make a little sanity check for presence of coordinates, if needed
    */
    if(nAlt>0)
        {
        if(A[0].m[mi][0].r[0].a[0].nalt<nAlt) /* if the first atom doesn't have enough alt coords */
            {mywhine("In get_assembly_PDB_ATOM_lines, requested coordinates do not appear to exist.");}
        } 
    /*
        Grab the ATOM lines for each molecule
    */
    FM[iFM]=get_molecule_PDB_ATOM_lines(A[0].m[mi],rinit,ainit,savei,savei,'n',raltname,this_x);
    if(FM[iFM].n==0){mywhine("no lines found in FM from get_assembly_PDB_ATOM_lines ");}
    ntot+=FM[iFM].n;
    /* 
      get new values for rinit and ainit if not using 'n'
         must use L[n-2] because the last one is always "TER"
    */
    if(ainit!=-1)
        {
        aptr=&FM[iFM].L[FM[iFM].n-2][6]; /* serial begins index 6 (column 7) */
         tmp=strdup(aptr); /* 5 chars */
         tmp[6]='\0';
         sscanf(tmp,"%d",&ainit);
         tmp=NULL;
         ainit++;
        rptr=&FM[iFM].L[FM[iFM].n-2][22]; /* resSeq begins index 22 (column 23) */
         tmp=strdup(rptr);
         tmp[5]='\0';
         sscanf(tmp,"%d",&rinit);
         tmp=NULL;
         rinit++;
        }
    iFM++;
    } /* close loop over number of molecules */
    if(ainit!=-1) { ainit=rinit=1; }
    this_x++;
    } /* close loop over number of models */
/*
    Set up the space for the whole assembly's collected lines 
*/
if(nModels>1) {ntot += nModels*2;} /* for MODEL and ENDMDL lines */
FA.n=ntot;
FA.L=(char**)calloc(FA.n,sizeof(char*));
Li=0; /* used here to be the total number of lines */
for(Mi=0;Mi<nModels;Mi++)
    {
    iFM=Mi*A[0].nm;
    if(nModels>1) 
        { 
        FA.L[Li]=(char*)calloc(16,sizeof(char));
        sprintf(FA.L[Li],"MODEL     %4d\n",Mi+1);
        Li++;
        }
for(mi=0;mi<A[0].nm;mi++)
    {
    for(ai=0;ai<FM[iFM].n;ai++)
        {
        FA.L[Li]=strdup(FM[iFM].L[ai]);
        Li++;
        }
    deallocateFileslurp(&FM[iFM]);
    iFM++;
    } 
    if(nModels>1) 
        { 
        FA.L[Li]=(char*)calloc(8,sizeof(char));
        sprintf(FA.L[Li],"ENDMDL\n");
        Li++;
        }
    } /* close loop over number of models */
free(FM);

return FA;
}

fileslurp get_ensemble_PDB_ATOM_lines(ensemble *E,char isource,int savei,char raltname,int xs)
{
int mi,ri,ai,ainit,rinit,ntot,Li,nModels,Mi,this_x,nAlt,iFM;
fileslurp FE,*FM; 
char *aptr,*rptr,*tmp;
/* allocate savei space if not -1 and if not already allocated */
if(savei!=-1)
    {
    for(mi=0;mi<E[0].nm;mi++)
        {
        for(ri=0;ri<E[0].m[mi].nr;ri++)
            {
            if(E[0].m[mi].r[ri].ni==0) 
                {
                E[0].m[mi].r[ri].ni=savei+1;
                E[0].m[mi].r[ri].i=(int*)calloc(savei+1,sizeof(int));
                }
            else if (E[0].m[mi].r[ri].ni<savei+1)
                {
                E[0].m[mi].r[ri].ni=savei+1;
                E[0].m[mi].r[ri].i=\
                   (int*)realloc(E[0].m[mi].r[ri].i,(savei+1)*sizeof(int));
                }
            for(ai=0;ai<E[0].m[mi].r[ri].na;ai++)
                {
                if(E[0].m[mi].r[ri].a[ai].ni==0) 
                    {
                    E[0].m[mi].r[ri].a[ai].ni=savei+1;
                    E[0].m[mi].r[ri].a[ai].i=(int*)calloc(savei+1,sizeof(int));
                    }
                else if (E[0].m[mi].r[ri].a[ai].ni<savei+1)
                    {
                    E[0].m[mi].r[ri].a[ai].ni=savei+1;
                    E[0].m[mi].r[ri].a[ai].i=\
                       (int*)realloc(E[0].m[mi].r[ri].a[ai].i,(savei+1)*sizeof(int));
                    }
        } } }
    }
/* 
    As a quick error check, see if the first residue's altname exists, if requested. 
*/
if((tolower(raltname)=='y')&&(E[0].m[0].r[0].altname==NULL))
    {
    printf("\nYou requested an assembly PDB print using alternate residue names.\n");
    printf("But the first one has not been set.  Not checking others.\n");
    printf("This might cause trouble. (in get_ensemble_ATOM_lines)\n");
    }

/* 
  Figure out how many slurps we will need based on the coords to be printed.
*/
if(xs==-3)  /* print main plus all alternates */  
  { 
  nAlt=E[0].m[0].r[0].a[0].nalt; /* for error-checking, grab number alts for first atom */
  nModels=nAlt+1; 
  this_x=-1;
  }
else if (xs==-2) /* print only all alternates */
  {
  nAlt=E[0].m[0].r[0].a[0].nalt; /* for error-checking, grab number alts for first atom */
  nModels=nAlt; 
  this_x=0; 
  }
else if (xs==-1) /* print only the main set of coordinates */
  {
  nAlt=0; 
  nModels=1; 
  this_x=-1; 
  }
else /* print only one alternate set */
  {
  nAlt=xs+1; 
  nModels=1; 
  this_x=xs; 
  }

FM=(fileslurp*)calloc((nModels*E[0].nm),sizeof(fileslurp));
if(tolower(isource)=='n') ainit=rinit=-1;
else ainit=rinit=1;

ntot=0; /* the total number of lines returned by get_molecule_PDB_ATOM_lines */

for(Mi=0;Mi<nModels;Mi++)
    {
    iFM=Mi*E[0].nm;
for(mi=0;mi<E[0].nm;mi++)
    {
    /*
        Make a little sanity check for presence of coordinates, if needed
    */
    if(nAlt>0)
        {
        if(E[0].m[mi].r[0].a[0].nalt<nAlt) /* if the first atom doesn't have enough alt coords */
            {mywhine("In get_assembly_PDB_ATOM_lines, requested coordinates do not appear to exist.");}
        } 
    /*
        Grab the ATOM lines for each molecule
    */
    FM[iFM]=get_molecule_PDB_ATOM_lines(&E[0].m[mi],rinit,ainit,savei,savei,'n',raltname,this_x);
    if(FM[iFM].n==0){mywhine("no lines found in FM from get_ensemble_PDB_ATOM_lines ");}
    ntot+=FM[iFM].n;
    /* 
      get new values for rinit and ainit if not using 'n'
         must use L[n-2] because the last one is always "TER"
    */
    if(ainit!=-1){
        aptr=&FM[iFM].L[FM[iFM].n-2][6]; /* serial begins index 6 (column 7) */
         tmp=strdup(aptr);
         tmp[6]='\0';
         sscanf(tmp,"%d",&ainit);
         tmp=NULL;
         ainit++;
        }
    if(rinit!=-1){
        rptr=&FM[iFM].L[FM[iFM].n-2][22]; /* resSeq begins index 22 (column 23) */
         tmp=strdup(rptr);
         tmp[5]='\0';
         sscanf(tmp,"%d",&rinit);
         tmp=NULL;
         rinit++;
        }
    iFM++;
    } /* close loop over number of molecules */
    this_x++;
    } /* close loop over number of models */
/*
    Set up the space for the whole assembly's collected lines 
*/
if(nModels>1) {ntot += nModels*2;} /* for MODEL and ENDMDL lines */
FE.n=ntot;
/*printf("ntot is %d\n",ntot);*/
FE.L=(char**)calloc(FE.n,sizeof(char*));
Li=0; /* used here to be the total number of lines */
for(Mi=0;Mi<nModels;Mi++)
    {
    iFM=Mi*E[0].nm;
    if(nModels>1) 
        { 
        FE.L[Li]=(char*)calloc(16,sizeof(char));
        sprintf(FE.L[Li],"MODEL     %4d\n",Mi+1);
        Li++;
        }
for(mi=0;mi<E[0].nm;mi++)
    {
    for(ai=0;ai<FM[iFM].n;ai++)
        {
        FE.L[Li]=strdup(FM[iFM].L[ai]);
        Li++;
        }
    deallocateFileslurp(&FM[iFM]);
    iFM++;
    } 
    if(nModels>1) 
        { 
        FE.L[Li]=(char*)calloc(8,sizeof(char));
        sprintf(FE.L[Li],"ENDMDL\n");
        Li++;
        }
    } /* close loop over number of models */
free(FM);

return FE;
}


const char *
get_PDB_line_for_ATOM(atom *a, residue *r, int ai, int ri, int asave, char raltname,int xs)
{
char *line, tmp[81];
int i,rhere=0,ahere=0;
if((xs>=0)&&(a[0].xa==NULL))
    {mywhine("In get_PDB_line_for_ATOM, set alternate coordinate, but a.xa is NULL.");}
line=(char*)calloc(82,sizeof(char)); /* 80 cols plus newline and termination */
/*  
pdb_a[2].b[1].c[0]	=	6;   RECORD NAME	(1-6)
*/
sprintf(line,"ATOM  ");
/*
pdb_a[2].b[1].c[1]	=	5;   serial		(7-11)

    This gets saved to a[0].i[asave] as well for LINK and CONECT cards later
    a[0].i[asave] space must be preallocated.
*/
if(ai==-1) ahere=a[0].n;
else ahere=ai;
if(ahere>99999){mywhine("Can not have atom serial > 99999 in pdb write utils.\n");}
strcat(line,get_float_string((double)ahere, 'r', 5, 0));
if(asave>=0){a[0].i[asave]=ahere;}
/*
pdb_a[2].b[1].c[2]	=	1;   N/A 		(12)
*/
strcat(line," ");
/*
pdb_a[2].b[1].c[3]	=	4;   name		(13-16)

    If there is an element (a.E) defined, use that to justify 
    else, pad initially with a space unless (a) the name starts with
    a number or (b) the name is already 4 chars.

    Ideally, figure out what element we have and set the beginning
    of it's designator at space 2, but... code that later...

    In the meantime, the calling function will have to set the name
    to a pdb-appropriate value if this function can't do it.
*/
tmp[0]='\0';
if((strlen(a[0].N)>=4)||(isdigit(a[0].N[0]!=0))){ strcpy(tmp,a[0].N); }
else 
    {
    if(a[0].N!=NULL){
	for(i=0;i<strlen(a[0].N);i++){ 
            if(a[0].N[i]==a[0].N[0]){break;} 
            }
	if(i<strlen(a[0].N)){ /* we have an element name here */
		if(i==0){sprintf(tmp," %s",a[0].N);}
                else{strcpy(tmp,a[0].N);}
		}
	}
    else {sprintf(tmp," %s",a[0].N);}
    }
    
strcat(line,get_char_string(tmp,'l',4));
/*
pdb_a[2].b[1].c[4]	=	1;   altLoc 		(17)

    This is currently not supported in the library (20110410)
*/
strcat(line," ");
/*
pdb_a[2].b[1].c[5]	=	3;   resName 		(18-20)
*/
if(tolower(raltname)=='y')
    {
    if(r[0].altname!=NULL) strcat(line,get_char_string(r[0].altname,'l',3));
    else strcat(line,"???");
    }
else{strcat(line,get_char_string(r[0].N,'l',3));}
/*
pdb_a[2].b[1].c[6]	=	1;   N/A 		(21)
*/
strcat(line," "); 
/*
pdb_a[2].b[1].c[7]	=	1;   chainID        	(22)
*/
if((a[0].cID!='\0')&&(a[0].cID!=NULL)){strncat(line,a[0].cID,1);}
else if((r[0].cID!='\0')&&(a[0].cID!=NULL)){strncat(line,r[0].cID,1);}
else{strcat(line," ");}
/*
pdb_a[2].b[1].c[8]	=	4;   resSeq		(23-26)
*/
if(ri==-1) rhere=r[0].n;
else rhere=ri;
if(rhere>9999){mywhine("Can not have residue serial > 9999 in pdb write utils.\n");}
strcat(line,get_float_string((double)rhere, 'r', 4, 0));
/*
pdb_a[2].b[1].c[9]	=	1;   iCode		(27)
*/
if(r[0].IC=='\0'){strcat(line," ");}
else{strncat(line,r[0].IC,1);}
/*
pdb_a[2].b[1].c[10]	=	3;   N/A		(28-30)
*/
strcat(line,"   ");
/*
pdb_a[2].b[1].c[11]	=	8;   x			(31-38)
*/
if(xs==-1) strcat(line,get_float_string((double)a[0].x.i, 'r', 8, 3));
else strcat(line,get_float_string((double)a[0].xa[xs].i, 'r', 8, 3));
/*
pdb_a[2].b[1].c[12]	=	8;   y			(39-46)
*/
if(xs==-1) strcat(line,get_float_string((double)a[0].x.j, 'r', 8, 3));
else strcat(line,get_float_string((double)a[0].xa[xs].j, 'r', 8, 3));
/*
pdb_a[2].b[1].c[13]	=	8;   z			(47-54)
*/
if(xs==-1) strcat(line,get_float_string((double)a[0].x.k, 'r', 8, 3));
else strcat(line,get_float_string((double)a[0].xa[xs].k, 'r', 8, 3));
/*
pdb_a[2].b[1].c[14]	=	6;   occupancy		(55-60)

    This currently doesn't have a computational equivalent. 
    Perhaps it could represent an rmsd at standard ambient 
    conditions or something.  But, for now, we just assign a
    value of 1.00 so that certain other programs will accept
    the file.

*/
strcat(line,"  1.00");
/*
pdb_a[2].b[1].c[15]	=	6;   tempFactor		(61-66)

    This currently doesn't have a computational equivalent. 
    We just set it to zero.
*/
strcat(line,"  0.00");
/*
pdb_a[2].b[1].c[16]	=	10;   N/A 		(67-76)
*/
strcat(line,"          ");
/*
pdb_a[2].b[1].c[17]	=	2;   element        
*/
/*printf("a[0].E is >>>%s<<<\n",a[0].E);*/
if(a[0].E==NULL){strcat(line,"  ");}
else
  {
  if(strlen(a[0].E)==1)
    {
    strcat(line," ");
    strcat(line,a[0].E);
    }
  else{strncat(line,a[0].E,2);}
  }

/*
pdb_a[2].b[1].c[18]	=	2;   charge         

    This has little meaning unless writing a "pdbq" (e.g., autodock)
    file, which breaks standard format.  We just leave blank.
*/
/*printf("the line (%d chars long) is >>>%s<<<\n",strlen(line),line);*/
strcat(line,"  \n");

if(strlen(line)>81){mywhine("something went wrong in get_PDB_line_for_ATOM");}


return (const char *)line;
}

void make_ATOM_HETATM(char *pdbline){
if(strlen(pdbline)<6){mywhine("pdbline too short in make_ATOM_HETATM");}
if((pdbline[0]=='A')&&(pdbline[1]=='T')&&(pdbline[2]=='O')&&(pdbline[3]=='M')){
  pdbline[0]='H';
  pdbline[1]='E';
  pdbline[2]='T';
  pdbline[3]='A';
  pdbline[4]='T';
  pdbline[5]='M';
  }
return;
}

/*
   ri and ainit are the numbers to use for the atom serial and 
   residue resSeq.  If set to -1, they will use the "n" values
   stored in the structure.  
*/
fileslurp 
get_residue_PDB_ATOM_lines(residue *r,int ri,int ainit,int rsave, int asave,char raltname, int xs)
{
int i,aihere,rihere;
fileslurp FR;
FR.n=r[0].na;
FR.L=(char**)calloc(FR.n,sizeof(char*));
if(ri==-1){rihere=r[0].n;}
else{rihere=ri;}
if(ainit>=0){aihere=ainit;}
/*
    r[0].i[rsave] and each a[i].i[asave] should be pre-allocated
	unless they are set to -1, and they won't be used
*/
if(rsave>=0){r[0].i[rsave]=rihere;}
for(i=0;i<FR.n;i++){
    if(ainit==-1){aihere=r[0].a[i].n;} 
    /* call get_PDB_line_for_ATOM for each atom */
    FR.L[i]=strdup(get_PDB_line_for_ATOM(&r[0].a[i], r, aihere, rihere, asave, raltname, xs));
    if(ainit>=0){aihere++;}
    }
return FR;
}

fileslurp 
get_molecule_PDB_ATOM_lines(molecule *mol,int rinit,int ainit,int rsave,int asave,
	char oneres, char raltname, int xs)
{
int i,j,b,aihere,rihere,*resconnects,ntot=0;
residue *this_r;
atom *this_a;
fileslurp *FR,FM;
char bondinfocomplete='y';

FR=(fileslurp*)calloc(mol[0].nr,sizeof(fileslurp));
resconnects=(int*)calloc(mol[0].nr,sizeof(int));

if(ainit==-1){aihere=-1;}
else{aihere=ainit;}
if(rinit>=0){rihere=rinit;}
if(tolower(oneres)=='y')
    {
    if(rinit==-1){rihere=mol[0].r[0].n;}
    else{rihere=rinit;}
    rinit=-2;
    }

for(i=0;i<mol[0].nr;i++){
    this_r=&mol[0].r[i];
    if(rinit==-1){rihere=this_r[0].n;} 
    FR[i]=get_residue_PDB_ATOM_lines(this_r,rihere,aihere,rsave,asave,raltname,xs);
    ntot+=FR[i].n;
    if(ainit>=0){aihere+=this_r[0].na;}
    for(j=0;j<this_r[0].na;j++){
        this_a=&this_r[0].a[j];
        if((this_r[0].na>1)&&(this_a[0].nmb==0)){bondinfocomplete='n';}
        for(b=0;b<this_a[0].nmb;b++){
            if(this_a[0].mb[b].s.r!=this_a[0].mb[b].t.r){
                resconnects[i]++;
                }
            }
        } 
    if(rinit>=0){rihere++;}
    }

/*
    Until the residue-level connection tree gets set, there is no foolproof
    or elegant way to do this next bit.  In fact, until the PDB folks adopt
    an official "end of branch" card, there isn't even a legitimate way to
    do it.  The following assumes the first residue is the connection tree
    start point -- that is, if there are two+ connections to other residues,
    then there is a branch.  All other residues are branch points if they
    have three or more connections to other residues.
*/
if(bondinfocomplete=='n')
    {
    printf("\nBonding info not complete or not present in write PDB utils.\n");
    printf("\tCannot set LINK and CONECT cards or TER cards within molecules:\n");
    }
/* find the number of TER cards */
b=0; /* the number of TER cards needed */
if(resconnects[0]>1){b++;}
for(i=1;i<mol[0].nr;i++)
    { if(resconnects[i]>2){ b++; } }
b++; /* for the final TER card */
FM.n=ntot+b;
FM.L=(char**)calloc(FM.n,sizeof(char*));
b=0; /* now, it's the total number of lines */
for(i=0;i<mol[0].nr;i++)
    {
    for(j=0;j<FR[i].n;j++)
        {
	FM.L[b]=strdup(FR[i].L[j]);        
        b++;
        }
    if(i==0)
        {
        if(resconnects[0]>1)
            {
            FM.L[b]=strdup("TER   \n");
            b++;
            }
        }
    else
        {
        if(resconnects[i]>2)
            {
            FM.L[b]=strdup("TER   \n");
            b++;
            }
        }
    }
if(b>(FM.n-1)){mywhine("something went wrong in get_molecule_PDB_ATOM_lines");}
FM.L[b]=strdup("TER   \n");

return FM;
}




void outputMolPDB(molecule* mol, char* filename)
{
 int i;
 FILE *file;
 fileslurp FM; 
 FM=get_molecule_PDB_ATOM_lines(mol,1,1,-1,-1,'n','n',-1);
 if(FM.n==0){mywhine("no lines found in FM from outputMolPDB");}
 file = myfopen(filename,"w");
 for(i=0;i<FM.n;i++)
  {
  if(FM.L[i]==NULL){mywhine("found null line in FM from outputMolPDB");}
  if(strncmp("TER",FM.L[i],3)!=0)
   {
   fprintf(file,"%s",FM.L[i]);
   }
  }
 fclose(file);
}


void outputAsmblPDB(assembly* asmbl, char* filename)
{
 int i,mi,rinit=0,ainit=0;
 FILE* file;
 char *tmp,*rptr,*aptr;
 fileslurp *FM; 
 FM=(fileslurp*)calloc(asmbl[0].nm,sizeof(fileslurp));
 rinit=ainit=1;
 for(mi=0;mi<asmbl[0].nm;mi++)
  {
  FM[mi]=get_molecule_PDB_ATOM_lines(asmbl[0].m[mi],rinit,ainit,-1,-1,'n','n',-1);
  if(FM[mi].n==0){mywhine("no lines found in FM from outputAsmblPDB");}
  /* 
    get new values for rinit and ainit 
       must use L[n-2] because the last one is always "TER"
  */
  aptr=&FM[mi].L[FM[mi].n-2][6]; /* serial begins index 6 (column 7) */
   tmp=strdup(aptr);
   tmp[6]='\0';
   sscanf(tmp,"%d",&ainit);
   ainit++;
   free(tmp);
  rptr=&FM[mi].L[FM[mi].n-2][22]; /* resSeq begins index 22 (column 23) */
   tmp=strdup(rptr);
   tmp[5]='\0';
   sscanf(tmp,"%d",&rinit);
   rinit++;
   free(tmp);
  }
 file = myfopen(filename,"w"); 
 for(mi=0;mi<asmbl[0].nm;mi++)
  {
  for(i=0;i<FM[mi].n;i++)
   {
   if(FM[mi].L[i]==NULL){mywhine("found null line in FM from outputAsmblPDB");}
   if(strncmp("TER   ",FM[mi].L[i],6)!=0)
    {
    fprintf(file,"%s",FM[mi].L[i]);
    }
   }
  fprintf(file,"TER   \n");
  deallocateFileslurp(&FM[mi]);
  }
 free(FM);
 fclose(file);
}

