/*Author: Michael Tessier with changes by many others 
	including BLFoley, possibly also Spandana Makeneni,
	Anita Nivedha, Matthew Tessier */
#include <load_pdb.h>
assembly* load_pdb(char* file_name)
{
  char* temp;
  WARN=1;SILENT=1;INWC=0;DATE=0;UNCUT=1;UNx=1;UNy=1;UNz=1;
  temp = suf; strcpy(temp,".pdb");
  temp = sufc; strcpy(temp,"_change.txt"); ACT='0';
  DEBUG=-1;LASTRES=-1;LASTOKX=0;LASTOKY=0;LASTOKZ=0;
  UNCTOL=0;CRYX=0;CRYY=0;CRYZ=0;LASTX=0;LASTY=0;LASTZ=0;
  int ma,ra;
  assembly* asmbl;
  IN = myfopen(file_name, "r");
  char find_LINK='n',find_CONECT='n';
  /*
  init_struct() scans pdb file for number of lines
  allocates memory for line structures initializes 
  structures that contain line formats

  This, among other PDB functions, needs to be made less "global".
  */
  init_struct();
  rewind(IN);
  /*Read in lines to ln structure */
  //printf("Reading in %s...\n",file_name);
  for(ma=0;ma<INWC;ma++){
  if(DEBUG>=1){printf("made it to here main line loop %d ...\n", ma);} 
	rwm_line(ma+1); 
	if(strncmp(ln[ma].f[0].c,"CONECT",6)==0) find_CONECT='y';
	if(strncmp(ln[ma].f[0].c,"LINK",4)==0) find_LINK='y';
	}

  /*Determine the number of molecules in the pdb */
  //printf("There are %d molecule(s)\n",howManyMolecules());
  asmbl = getAssembly();
  set_assembly_molindexes(asmbl);
  if(find_CONECT=='y')
    { 
    set_assembly_atom_molbonds_from_PDB_CONECT(asmbl, ln, file_name, INWC); 

/* 
 * dropping this because apparently not working...
    for(ma=0;ma<asmbl[0].nm;ma++)
        {
        for(ra=0;ra<asmbl[0].m[ma][0].nr;ra++)
            {
            set_residue_atom_nodes_from_bonds(&asmbl[0].m[ma][0].r[ra]);
            }
        set_molecule_atom_nodes_from_bonds(asmbl[0].m[ma]);
*/
        /*
        The atom-level CONECT cards are more trusted than LINK cards,
        so they are used unless there are no CONECT cards.
        */
/*
 * also this
        set_molecule_residue_molbonds(asmbl[0].m[ma]);
        set_molecule_residue_nodes_from_bonds(asmbl[0].m[ma]);
        }
*/
    }
/* and this isn't even written
  if((find_LINK=='y')&&(find_CONECT=='n'))
    { 
    set_assembly_residue_molbonds_from_PDB_LINK(asmbl, ln, file_name, INWC); 
    }
*/
  //printf("PDB Information Successfully Read.\n");
  free(ln);
  return asmbl;
}

molecule load_pdb_from_slurp(fileslurp in){
 int i,hldr;
 char* temp;
 init_struct_slurp(in);
 for(i = 0; i < in.n; i++){	//Assign data to the rwm line
  rwm_line_char(in.L[i],i+1);
  if(isAtom((ln+i))){		//And determine the smallest residue
   temp  = (*(ln+i)).f[8].c; sscanf(temp,"%d",&hldr);
  }
 }
 INWC = in.n;
/*
printf("INWC = %d\n",INWC);
*/
 return (*getMolecule());
}

void set_assembly_atom_molbonds_from_PDB_CONECT(assembly *A, linedef *ln, const char *file_name, int nlines)
{
int ai,bi,li,ncurrent=-1,this_n,that_n,list_loc,na=A[0].na;
char badorder='n';
molindex_set found_a;
int this_ai,alist[na],bonds[na],blist[na],llist[na];
int conect_first=-1,conect_last=-1,tmpi,tmpj;
int checkNum=34685;
/*
    0.  Check that atoms are all in order.  Complain if not, but try anyway.
*/
/*
printf("A[0].na is %d\n",A[0].na);
printf("0. A[0].a[checkNum][0].n is %d and (*A).a[checkNum+1][0].n is %d\n",A[0].a[checkNum][0].n,A[0].a[checkNum+1][0].n);
*/
 for(ai=0;ai<A[0].na;ai++){
  if(A[0].a[ai][0].n<=ncurrent) 
	{ 

printf("A[0].a[%d][0].n (%d) is less than ncurrent (%d), and this concerns us\n",ai,A[0].a[ai][0].n,ncurrent);

	badorder = 'y'; 
	}
/*
printf("ncurrent is %d ; n is %d \n",ncurrent,A[0].a[ai][0].n); 
if((A[0].a[ai][0].n-ncurrent)!=1) printf("ncurrent (%d) differs from n (%d) by %d\n",ncurrent,A[0].a[ai][0].n,A[0].a[ai][0].n-ncurrent);
*/
	ncurrent=A[0].a[ai][0].n;
	}
/*
printf("1. A[0].a[checkNum][0].n is %d and (*A).a[checkNum+1][0].n is %d\n",A[0].a[checkNum][0].n,A[0].a[checkNum+1][0].n);
*/
 if(badorder=='y')
  {
  printf("\nAtoms serials (numbers) in file %s are not in increasing order.\n",file_name);
  printf("This might cause bonding problems.\n");
  }
for(ai=0;ai<A[0].na;ai++) 
  { 
  alist[ai]=-1; 
  blist[ai]=-1; 
  llist[ai]=-1; 
  bonds[ai]=0; 
  }
/*
    1.  Find numbers for atoms and number of bonds for each.
*/ 
/*
printf("2. A[0].a[checkNum][0].n is %d and (*A).a[checkNum+1][0].n is %d\n",A[0].a[checkNum][0].n,A[0].a[checkNum+1][0].n);
*/
ncurrent=0;
 for(li=0;li<nlines;li++)
  {
  if(strcmp(ln[li].f[0].c,"CONECT")==0)
   { 
   if(conect_first==-1) conect_first=li;
   sscanf(ln[li].f[1].c,"%d",&this_n);
   list_loc=-1;
   for(ai=0;ai<ncurrent;ai++) { if(alist[ai]==this_n) list_loc=ai; }
//printf("this_n is %d ; list_loc is %d \n",this_n,list_loc); 
   if(list_loc<0)
    { 
    list_loc=ncurrent;
    ncurrent++;
    alist[list_loc]=this_n;
    llist[list_loc]=li;
    }
//printf(" ---  list_loc is %d \n",list_loc); 
   for(bi=2;bi<15;bi++)
    {
    if(ln[li].f[bi].c==NULL) break;
    tmpi=sscanf(ln[li].f[bi].c,"%d",&tmpj); 
    if(tmpi<=0) { continue; }
//printf("\t\ttmpi is %d ; tmpj is %d\n",tmpi,tmpj); 
    bonds[list_loc]++;
    }
//printf("\t\tbonds[%d] is  %d\n",list_loc,bonds[list_loc]); 
    conect_last=li;
   }
  }
/*
printf("3. A[0].a[checkNum][0].n is %d and (*A).a[checkNum+1][0].n is %d\n",A[0].a[checkNum][0].n,A[0].a[checkNum+1][0].n);
*/
/*printf("conect first/last are:  %d   %d  \n",conect_first,conect_last); */
/*
    2.  Set atom-level connections.
*/ 
  for(ai=0;ai<na;ai++) 
   {
   if(alist[ai]==-1) continue;
   for(ncurrent=0;ncurrent<na;ncurrent++) { blist[ncurrent]=-1; }
   found_a=find_assembly_top_level_atoms_by_n(A,alist[ai]);
/*printf("found_a.nP is %d and P[0] i-m-r-a is %d-%d-%d-%d n=%d N=%s\n",found_a.nP,\
	found_a.P[0].i,found_a.P[0].m,found_a.P[0].r,found_a.P[0].a,\
	A[0].m[found_a.P[0].m][0].r[found_a.P[0].r].a[found_a.P[0].a].n,A[0].m[found_a.P[0].m][0].r[found_a.P[0].r].a[found_a.P[0].a].N);
*/
   if(found_a.nP==0)
    {
    printf("\nBonding specified in file %s for non-existent atom serial number %d\n",file_name,this_n);
    printf("Skipping this one and hoping for better next time.\n");
    continue;
    }
   if(found_a.nP>1)
    {
    printf("\nBonding specified in file %s for duplicated atom serial number %d\n",file_name,this_n);
    printf("Choosing first in list and moving on.\n");
    }
   this_ai=found_a.P[0].i;
   free(found_a.P);
   A[0].a[this_ai][0].nmb=bonds[ai];
   A[0].a[this_ai][0].mb=(molbond*)calloc(bonds[ai],sizeof(molbond));
   list_loc=0;
   for(li=conect_first;li<=conect_last;li++)
    {
    tmpi=sscanf(ln[li].f[1].c,"%d",&tmpj);
    if(tmpi<=0) {mywhine("trouble reading ln[li].f[1].c in set_assembly_atom_molbonds_from_PDB_CONECT");} 
    if(tmpj==alist[ai])
     {
     for(bi=2;bi<15;bi++)
      {
      if(ln[li].f[bi].c==NULL) break;
      tmpi=sscanf(ln[li].f[bi].c,"%d",&that_n); 
      if(tmpi<=0) { continue; }
      /* If there is an entry in this location */
      found_a=find_assembly_top_level_atoms_by_n(A,that_n);
      if(found_a.nP==0)
       {
       printf("\nBonding specified in file %s for non-existent atom serial number %d\n",file_name,that_n);
       printf("Skipping this one and hoping for better next time.\n");
       continue;
       }
      if(found_a.nP>1)
       {
        printf("\nBonding specified in file %s for duplicated atom serial number %d\n",file_name,that_n);
        printf("Choosing first in list and moving on.\n");
        } 
       A[0].a[this_ai][0].mb[list_loc].s=A[0].a[this_ai][0].moli;
       A[0].a[this_ai][0].mb[list_loc].t=found_a.P[0];
       free(found_a.P);
       list_loc++;
       if(list_loc>A[0].a[this_ai][0].nmb){mywhine("list_loc>A[0].a[this_ai][0].nmb in set_assembly_atom_molbonds_from_PDB_CONECT");}
       }
     }
    }
   if(list_loc<A[0].a[this_ai][0].nmb){mywhine("list_loc<A[0].a[this_ai][0].nmb in set_assembly_atom_molbonds_from_PDB_CONECT");}
   }
/*
printf("4. A[0].a[checkNum][0].n is %d and (*A).a[checkNum+1][0].n is %d\n",A[0].a[checkNum][0].n,A[0].a[checkNum+1][0].n);
*/
return;
}

void set_assembly_residue_molbonds_from_PDB_LINK(assembly *asmbl, linedef *ln, const char *file_name, int nlines)
{
mywhine("The function set_assembly_residue_molbonds_from_PDB_LINK has not been written yet.\n");
}


fileslurp isolateInputPDB(fileslurp S){
 int i,itr;
 char* plc;
 fileslurp O; O.n = 0; itr = 0;
 for(i = 0; i < S.n; i++){
  if(strstr(S.L[i],"INPUT-PDBQ: ATOM")!=NULL ||		//These cover the two types of
     strstr(S.L[i],"INPUT-PDBQ: HETATM")!=NULL ||	//AD3 style input data
     strstr(S.L[i],"INPUT-LIGAND-PDBQT: ATOM")!=NULL ||	//These cover the two types of
     strstr(S.L[i],"INPUT-LIGAND-PDBQT: HETATM")!=NULL)	//AD4 style input data
   O.n++;
 }
 O.L = (char**)calloc(O.n,sizeof(char*));
 for(i = 0; i < S.n; i++){
  if(strstr(S.L[i],"INPUT-PDBQ: ATOM")!=NULL ||
     strstr(S.L[i],"INPUT-PDBQ: HETATM")!=NULL ||
     strstr(S.L[i],"INPUT-LIGAND-PDBQT: ATOM")!=NULL ||
     strstr(S.L[i],"INPUT-LIGAND-PDBQT: HETATM")!=NULL){
   if(strstr(S.L[i],"INPUT-PDBQ: ")!=NULL)
    plc = S.L[i]+12;	//Based on what starts the line, duplicate the line minus
   else			//The first couple of words
    plc = S.L[i]+20;
   O.L[itr] = strdup(plc); itr++;
  }
 }
 for(i = 0; i < S.n; i++)
  free(S.L[i]);
 free(S.L);
 return O;
}

fileslurp isolateDockedPDB(fileslurp S){
 int i,itr;
 char* plc;
 fileslurp O; O.n = 0; itr = 0;
 for(i = 0; i < S.n; i++){
  if(strstr(S.L[i],"DOCKED: ATOM")!=NULL)
   O.n++;
 }
 O.L = (char**)calloc(O.n,sizeof(char*));
 for(i = 0; i < S.n; i++){
  if(strstr(S.L[i],"DOCKED: ATOM")!=NULL){
   plc = S.L[i]+8;
   O.L[itr] = strdup(plc); itr++;
  }
 }
 for(i = 0; i < S.n; i++)
  free(S.L[i]);
 free(S.L);
 return O;
}


int howManyMolecules()
{
 int i; int counter = 0;
 char in_molecule_switch = 'n';
 int check_status=0;
/* Adding a little robustness to this function.
	20100727, BLFoley.  
   In particular, making a check for the possibility that there
	are no molecules found.  Also fixing issue due to there
	being multiple TER, CONECT, etc., cards in a row. */ 
 for(i = 0; i < INWC; i++) {
	check_status = isAtom((ln+i)); /* is this an atom? (ATOM or HETATM) */
	if(in_molecule_switch == 'n' && check_status == 1){
		in_molecule_switch = 'y'; /* we have entered a molecule */
		counter++;
		}
	check_status = endOfMol((ln+i)); /* If a card that usually ends a 
		molecule: TER, LINK, CONECT, END, ENDMDL, MASTER */
	if(check_status == 1){ in_molecule_switch = 'n'; }
	} 
 return counter;
}

assembly* getAssembly()
{
  int molNum = howManyMolecules();
  int i = 0;int j = 0;int k = 0;int ka = 0; int la=0;int atom_num;
  int init=0, next_mol_line=0, ri=0,natoms=0,resNum,Kref=0;
  char in_molecule_switch = 'n';
  char is_same_residue = 'y';
  double x,y,z;
  char atmName[5], atmElem[3], resName[4],IC=' ',cID=' ',*temp; 
  assembly* asmbl = (assembly*) calloc (1 , sizeof(assembly));
  (*asmbl).nm = molNum;
  molecule *curMol=NULL; residue* curRes=NULL; 
  atom* curAtm;// atom* atm;
  /* made further changes for molecule double pointers in assembly
	starting on 20100723.  BLFoley */ 

  (*asmbl).m = (molecule**)calloc(molNum,sizeof(molecule*));//added by MNT on 20080806
  for(init=0;init<molNum;init++){ 
	(*asmbl).m[init] = (molecule*) calloc (1,sizeof(molecule)); 
	(*asmbl).m[init][0].nr=-1; /* initializing as empty/not-seen. */ 
	}
  (*asmbl).nr = 0; // initializing 
  (*asmbl).r = (residue**) calloc(1, sizeof(residue*));
  (*asmbl).na = 0; // initializing 
  (*asmbl).a = (atom**) calloc(1, sizeof(atom*));
  i=0; 
  strcpy(atmName,"");
  strcpy(resName,"");
  strcpy(atmElem,"");
/*printf("First allocate of residues and atoms:\n");*/
/*printf("\t molecule %d, init is %d; residue %d -- nm is %d, A.nr is %d, A.na is %d\n",j,init,k,(*asmbl).nm,(*asmbl).nr,(*asmbl).na);*/
/* 
	i = line number
	j = molecule number
	k = residue number (within the molecule, not absolute)
	ka = absolute residue number 
        (*curRes).ni = current atom to place (set back to zero when done)
	la = absolute atom number 
*/

in_molecule_switch = 'n';
for(i = 0; i < INWC; i++) {
	if(j==molNum){break;}
	if(in_molecule_switch == 'n'){ /* If we are not in a molecule. */
		while( (i<(INWC)) && (isAtom(ln+i)==0) ){ i++; } /* advance to first ATOM/HETATM line */
		if(i==INWC){break;}
		curMol = &asmbl[0].m[j][0];
		(*curMol).Des=(char*)calloc(2,sizeof(char));
		(*curMol).Des[0]=(*(ln+i)).f[7].c[0]; (*curMol).Des[1]='\0';
		(*curMol).nr = findTotalResidue(i);
		(*curMol).r = (residue*) calloc ((*curMol).nr,sizeof(residue));
		/*Initialize the residue numbers to -1 so the program knows they're not used yet. */
		for(ri=0;ri<(*curMol).nr;ri++){
			initialize_residue(&curMol[0].r[ri]);
			(*curMol).r[ri].n = -1;
			(*curMol).r[ri].ni = 0;
			(*curMol).r[ri].moli.m=j; 
			(*curMol).r[ri].moli.r=ri; 
			(*curMol).r[ri].moli.a=-1;
			}
		/*Get the total atoms in the molecule while we're here. */
		next_mol_line=getResInfo((*curMol).r,i);
		natoms=0;
		for(ri = 0; ri < (*curMol).nr; ri++){
			(*curMol).r[ri].a = (atom*) calloc ((*curMol).r[ri].na,sizeof(atom));
			natoms+=(*curMol).r[ri].na;
/*
printf("na is %d for residue %d(ka=%d)\n",(*curMol).r[ri].na,ri,ka);
*/
			}
		/* Allocate the top-level atoms and residues for this molecule */
		(*asmbl).nr+=(*curMol).nr;
		(*asmbl).r = (residue**) realloc ((*asmbl).r, (*asmbl).nr*sizeof(residue*));
		(*asmbl).na+=natoms; 
		(*asmbl).a = (atom**) realloc ((*asmbl).a, (*asmbl).na*sizeof(atom*));
		/* set the top-level residues */
		for(ri = 0; ri < (*curMol).nr; ri++) { (*asmbl).r[ka+ri] = &curMol[0].r[ri]; }
		in_molecule_switch = 'y';/* At this point, we should be at an ATOM/HETATM entry. */
		curRes = curMol[0].r;
		curAtm = (*curRes).a;
		ka+=curMol[0].nr;
		} 
	if(i==INWC){ /* If we ran out of lines... */
		if((*curMol).nr==-1){ /* If this is molecule does not currently contain residues */
			(*curMol).nr = 0;
			(*curMol).r = NULL;
			(*curMol).N = strdup("EMPTY_MOLECULE");
			if(j!=molNum){printf("Found %d molecules, but ran out of atoms before molecule %d.\n",molNum,j+1);
				printf("This is probably a very, very bad thing, but we're ignoring it for now.\n");}
			} 
		else if ((*curMol).nr>0){ /* if this molecule already contains residues */
			printf("Got to i=INWC while in a molecule.  Go fix code.\n"); }
		else { printf("Got to i=INWC for an uncoded value of nr.  Go fix code.\n");}
		}
	while( (i<INWC) && (in_molecule_switch=='y') ) 
		{ /* while we are in a molecule */
		/* 2014-08-29 BLFoley -- fixed logic issues that affected single-atom residues */
		if(endOfMol((ln+i)) == 1){
			in_molecule_switch='n';
			if(next_mol_line!=i) {printf("i is %d and next_mol_line is %d ... should match.\n",i,next_mol_line);}
			break;
			}
		else {
			if(isAtom((ln+i)) == 1) {
/*
printf("\nTOP:  k = %d ; ka = %d ; j = %d ; i = %d\n",k,ka,j,i);
*/
				/* Figure out which residue the atom goes in */
				is_same_residue='y';
				temp  = ln[i].f[8].c; sscanf(temp,"%d",&resNum);
				temp = ln[i].f[5].c; sscanf(temp,"%s",resName);
				IC = ln[i].f[9].c[0];
				cID = ln[i].f[7].c[0];
/*
printf("0. resnum is %d ;    resname is >>%s<<     cID is >>%c<<           IC is >>%c<< \n",resNum,resName,cID,IC);
printf("   *curres.n is %d ; (*curRes).N is >>%s<< (*curRes).cID is >>%s<< (*curRes).IC is >>%s<< \n",(*curRes).n,(*curRes).N,(*curRes).cID,(*curRes).IC);
*/
				if((*curRes).n!=resNum){ is_same_residue='n';}
				if(strcmp((*curRes).N,resName)!=0){ is_same_residue='n';}
				if((*curRes).IC[0]!=IC){ is_same_residue='n';}
				if((*curRes).cID[0]!=cID){ is_same_residue='n';}
				if(is_same_residue=='n'){ /* try the rest of the residues */
				Kref=k+1;
/*
printf("checking other residues: curRes.na=%d ; i=%d ; Kref=%d \n",curRes[0].na,i,Kref);
*/
				for(k=Kref;k<(*curMol).nr;k++) {
					is_same_residue='y';
					curRes = &curMol[0].r[k];
/*
printf("checking the rest of the residues: curMol.nr=%d ; k=%d ; Kref=%d \n",curMol[0].nr,k,Kref);
printf("2. resnum is %d ;    resname is >>%s<<     cID is >>%c<<           IC is >>%c<< \n",resNum,resName,cID,IC);
printf("   *curres.n is %d ; (*curRes).N is >>%s<< (*curRes).cID is >>%s<< (*curRes).IC is >>%s<< \n",(*curRes).n,(*curRes).N,(*curRes).cID,(*curRes).IC);
*/
					if((*curRes).n!=resNum){ is_same_residue='n';}
					if(strcmp((*curRes).N,resName)!=0){ is_same_residue='n';}
					if((*curRes).IC[0]!=IC){ is_same_residue='n';}
					if((*curRes).cID[0]!=cID){ is_same_residue='n';}
/*
printf("\t2. Is same residue (N=%s  n=%d IC=%c cID=%c) = %c\n",(*curRes).N,(*curRes).n,(*curRes).IC[0],(*curRes).cID[0],is_same_residue);
*/
					if(is_same_residue=='y') { break; }
					}}


				if(is_same_residue=='n'){ /* try the first residues */
				for(k=0;k<(Kref-1);k++) {
					is_same_residue='y';
					curRes = &(*curMol).r[k];
/*
printf("checking the first residues: curMol.nr=%d ; k=%d ; Kref=%d \n",curMol[0].nr,k,Kref);
*/
					if((*curRes).n!=resNum){ is_same_residue='n';}
					if(strcmp((*curRes).N,resName)!=0){ is_same_residue='n';}
					if((*curRes).IC[0]!=IC){ is_same_residue='n';}
					if((*curRes).cID[0]!=cID){ is_same_residue='n';}

/*
printf("\t3. Is same residue (N=%s  n=%d IC=%c cID=%c) = %c\n",(*curRes).N,(*curRes).n,(*curRes).IC[0],(*curRes).cID[0],is_same_residue);
*/

					if(is_same_residue=='y') { break; }
					}}
				if(is_same_residue=='n'){mywhine("Can't find residue in load_pdb's getAssembly.");}
/*
printf("ni is %d and na is %d \n",(*curRes).ni,(*curRes).na);
*/
				if((*curRes).ni==(*curRes).na){mywhine("(*curRes).ni==(*curRes).na in getAssembly\n");}
				curAtm = &(*curRes).a[(*curRes).ni];
				asmbl[0].a[la] = &(*curRes).a[(*curRes).ni];
				if(curAtm[0].n!=0)
					{
					mywhine("curAtm[0].n!=0 in a place where it really should be.\n");
/*
printf("\nWe have a problem:\n");
printf("\tla=%d ; curRes[0].ni=%d ; k=%d ; ka=%d ; j=%d ; i=%d\n",la,curRes[0].ni,k,ka,j,i);
printf("ln[i].f[1].c=>>%s<<\n",ln[i].f[1].c);
printf("\tcurAtm[0].n=%d ; curRes[0].na=%d \n",curAtm[0].n,curRes[0].na);
printf("\tcurRes[0].n=%d ; curRes[0].IC=>>%s<< ; curRes[0].cID=>>%s<< \n",curRes[0].n,curRes[0].IC,curRes[0].cID);
fflush(stdout);
*/
					}

				//Getting the atom # //
				sscanf( (*(ln+i)).f[1].c ,"%d",&atom_num);
				//Getting the atom name //temp  = (*(ln+i)).f[3].c;
				if(((*(ln+i)).f[3].c==NULL)||((*(ln+i)).f[3].c[0]=='\0')){strcpy(atmName,"    ");}
				else{sscanf( (*(ln+i)).f[3].c ,"%s",atmName);}
				if(((*(ln+i)).f[17].c==NULL)||((*(ln+i)).f[17].c[0]=='\0')){strcpy(atmElem,"  ");}
				else{sscanf( (*(ln+i)).f[17].c ,"%s",atmElem); }
/*
printf("\t atmName is >>>%s<<< atmElem is >>>%s<<<\n",atmName,atmElem);
*/ 
				//Getting the X, Y and Z coordinates
				sscanf( (*(ln+i)).f[11].c ,"%lf",&x);
				sscanf( (*(ln+i)).f[12].c ,"%lf",&y);
				sscanf( (*(ln+i)).f[13].c ,"%lf",&z);   
				(*curAtm).n = atom_num;//set the atom #
				(*curAtm).N = strdup(atmName);//set the atom name
				(*curAtm).E = strdup(atmElem);//set the atom element
				// Record the chain identifier at the atom level
				if((*(ln+i)).f[7].c!=NULL){
					(*curAtm).cID=(char*)calloc(2,sizeof(char));
					(*curAtm).cID[0] = (*(ln+i)).f[7].c[0];
					(*curAtm).cID[1] = '\0';
					}
				else{strcpy((*curAtm).cID," ");}
/*
printf("\t the (*curAtm).cID is >>>%s<<< E is >>>%s<<<  the atom serial is %d\n",(*curAtm).cID,(*curAtm).E,(*curAtm).n); 
*/
				//set the atom's coordinates
				(*curAtm).x.i = x; (*curAtm).x.j = y; (*curAtm).x.k = z;
				(*curAtm).moli.m=j; (*curAtm).moli.r=k; (*curAtm).moli.a=(*curRes).ni;
				(*curRes).ni++;
				if(((*curRes).na==1)&&(k<(*curMol).nr)) { curAtm = (*curRes).a; }
				if((*curRes).ni==(*curRes).na) (*curRes).ni=0;
				la++;
				}
			}
		i++;
		} /* close while we are in a molecule */
	j++;
	}
/*
printf("ka=%d ; (*asmbl).nr=%d ; la=%d ; (*asmbl).na=%d\n",ka,(*asmbl).nr,la,(*asmbl).na);
*/
if((*asmbl).na != la){mywhine("(*asmbl).na != la in load_pdb's getAssembly");}
if((*asmbl).nr != ka){mywhine("(*asmbl).nr != ka in load_pdb's getAssembly");}

return asmbl;
}

int findTotalResidue(int start)
{
  int count=0;int i = start;int thisRes;//int curRes = 0;
  int inList,j;
  char *temp, icode=' ',chainid=' ';
  int  *known = (int*)calloc(1,sizeof(int));
  char *knownI = (char*)calloc(1,sizeof(char));
  char *knownID = (char*)calloc(1,sizeof(char));

//printf(" i is %d, and INWC is %d and ln[i].a-b are %d-%d\n",i,INWC,ln[i].a,ln[i].b);
//printf("isAtomln+i is %d\n",isAtom(ln+i));
//printf("endof Mol is %d\n",endOfMol(ln+i));  
while( (i!=INWC) && (isAtom((ln+i))==0) ) { i++; }

  while( (i != INWC) && (endOfMol((ln+i)) == 0) )
  {
    if(isAtom((ln+i)) == 1)
    {
     temp  = (*(ln+i)).f[8].c;
     icode = ln[i].f[9].c[0];  // we know this is only one character
     chainid = ln[i].f[7].c[0];  // we know this is only one character
     sscanf(temp,"%d",&thisRes);
     inList = 0;
     //The New Way
     for(j = 0; j < count; j++){
//printf("thisRes=%d icode=%c chainid=%c\n",thisRes,icode,chainid);
//printf("known=%d knownI=%c knownID=%c  j=%d\n",known[j],knownI[j],knownID[j],j);
      if((known[j]==thisRes)&&(knownI[j]==icode)&&(knownID[j]==chainid)){inList=1; break;}
//printf("still here\n");
     }
     if(!inList){
      count++;
      known = (int*)realloc(known,count*sizeof(int));
      known[count-1] = thisRes;
      knownI = (char*)realloc(knownI,count*sizeof(char));
      knownI[count-1] = icode;
      knownID = (char*)realloc(knownID,count*sizeof(char));
      knownID[count-1] = chainid;
//printf("known=%d knownI=%c knownID=%c  count-1=%d\n",known[count-1],knownI[count-1],knownID[count-1],count-1);
     }
    }
   i++;
  }
free(known);
free(knownI);
free(knownID);
//printf("count is %d\n",count);
  return count;
}

int endOfMol(linedef* line)
{
/* Originally, this checked for these cards:
	CONECT, LINK, TER.  
   Adding ENDMDL, MASTER and END.  20100727, BLFoley */
 switch((*line).a){ 
/*  if(((*line).a == 4 && (*line).b == 2)||
     ((*line).a == 2 && ((*line).b == 9||(*line).b == 3))) */
	case 0:
		switch ((*line).b) {
			case 1: /* END */
			case 3: /* MASTER */
				return 1;
			default: break;
			}
	case 2:
		switch ((*line).b) {
			case 3:  /* CONECT */
			case 9: /* LINK */
				return 1;
			default: break;
			}
	case 4:
		switch ((*line).b) {
			case 0: /* ENDMDL */
			case 2: /* TER */
				return 1;
			default: break; 
			}
	default: return 0;
	}
}

int isAtom(linedef* line)
{
  if(((*line).a == 2 || (*line).a == 3) && (*line).b == 1)
    return 1;
  else
    return 0;
}
int getResInfo(residue* res, int start)
{
int resNum,j;int i = start;
int count = 0, current_status=0;
residue* curRes = (res+0);
char name [10];char* temp;char* resName = name;
char IC=' ', cID=' ',is_same='y';

//printf("\n**1.  i is %d; atom name is %s; atom number is %s\n",i,(*(ln+i)).f[3].c,(*(ln+i)).f[1].c);

while( (current_status == 0) && (i<INWC) ){ 
	current_status = isAtom(ln+i); /* find out if we have an atom line  */
	i++; } /* if not, check the next line */
/*
printf("i is %d  and INWC is %d \n",i,INWC);
printf("current ststus is %d\n",current_status);
*/
if(i==INWC){ /* if we got to the end with no atoms */
//if current status is that we have an atom and i==1 then break
	if(current_status==1 && i==1){
	
	}
	else{
	(*curRes).N = strdup("EMPTY_RESIDUE");
	(*curRes).n = -100;
	(*curRes).na = 0;
	return INWC;
	}
	} 
i--; /* still here? decrement the counter to undo the while loop above */

while( (i != INWC)  && (endOfMol((ln+i)) == 0) ){ 
	if(isAtom(ln+i) == 1){
		temp  = ln[i].f[8].c; sscanf(temp,"%d",&resNum);
		temp = ln[i].f[5].c; sscanf(temp,"%s",resName);
		IC = ln[i].f[9].c[0];
		cID = ln[i].f[7].c[0];
		is_same='y';
		if(resNum != (*curRes).n) is_same='n';
		if((curRes[0].IC==NULL)||(IC!=curRes[0].IC[0])) is_same='n';
		if((curRes[0].cID==NULL)||(cID!=curRes[0].cID[0])) is_same='n';
/*
printf("resNum is %d and (*curRes).n is %d ; resName is >%s<<\n",resNum,(*curRes).n,resName);
printf("**2.  i is %d; atom name is %s; atom number is %s\n",i,(*(ln+i)).f[3].c,(*(ln+i)).f[1].c);
*/
		//If this is a different residue than the one in the previous line
		/*if((resNum != (*curRes).n)||(curRes[0].IC==NULL)||(IC!=curRes[0].IC[0]))*/
		if(is_same=='n')
			{
/*
printf("This residue is different from the previous.  Checking for broken residue.\n");
*/
			curRes = NULL;
			for(j = 0; j < count; j++){//Cycle through all of the found residues
				if(((*(res+j)).n==resNum)&&((*(res+j)).IC[0]==IC)&&((*(res+j)).cID[0]==cID)){ /* if already seen... */
					curRes = (res+j); break;} /* reset residue pointer */
				} 
			/* If it is not a known residue, assign it a new slot */
			if(curRes == NULL || curRes == 0x0){ curRes = (res+count); count++; }
			}
		if((*curRes).n < 0) { /* If this residue has not been found yet */
/*
printf("This residue has not been found yet.  Assigning values\n");
*/
			if(resName!=NULL){(*curRes).N = strdup(resName);}/* Set the residue name */
			else{(*curRes).N = strdup("   ");}
			curRes[0].IC=(char*)calloc(2, sizeof(char));
			(*curRes).IC[0] = IC;	/* Set the IC */
			(*curRes).IC[1] = '\0';	/* Set the IC */
			curRes[0].cID=(char*)calloc(2, sizeof(char));
			(*curRes).cID[0] = cID;	/* Set the cID */
			(*curRes).cID[1] = '\0';	/* Set the cID */
			(*curRes).n = resNum;		/* Set the actual residue number */
			(*curRes).na = 0;			/* ...and make sure the total # of atoms is 0 */
			} /* Otherwise we can assume all of these have already been set */
			
/*
printf("  -->  resNum is %d and (*curRes).n is %d, (*curRes).na is %d \n",resNum,(*curRes).n,(*curRes).na);
printf("**3.  i is %d; atom name is %s; atom number is %s\n",i,(*(ln+i)).f[3].c,(*(ln+i)).f[1].c);
*/
		(*curRes).na++;
/*
printf(" --> (curRes). n=%d, na=%d, N=%s\n",(*curRes).n,(*curRes).na,(*curRes).N);
printf("**4.  i is %d; atom name is %s; atom number is %s\n",i,(*(ln+i)).f[3].c,(*(ln+i)).f[1].c);
*/
		}//End if an atom 
	i++; //..and then incriment to the next line
	} //End loop through file
return i;
}

molecule* getMolecule(void){
  int molNum = 1;
  int i = 0;//int j = 0;int k = 0;
  int ntX,j,rI,next_mol_line;
  double dblY;
  char* temp;char name[6] = ""; //char* atmName = name;
  char icode=' ',chainid=' ';
  molecule* mol = (molecule*)calloc(molNum,sizeof(molecule));
  residue* curRes;// residue* res;
  atom* curAtm;// atom* atm;
  //Get all the information for the first molecule
  (*mol).nr = findTotalResidue(i);
  //printf("the number of residues are %d\n",(*mol).nr);
  (*mol).r = (residue*)calloc((*mol).nr,sizeof(residue));
  //Initialize the residue numbers to -1 so the program knows they're not used yet
  for(j = 0; j < (*mol).nr; j++){(*mol).r[j].n = -1;}
  next_mol_line=getResInfo((*mol).r,i);
//printf("next_mol_line is %d\n",next_mol_line);
  int l[(*mol).nr];
  for(i = 0; i < (*mol).nr; i++){
   	//printf("the number of atoms are %d\n",(*mol).r[i].na);
	(*mol).r[i].a = (atom*)calloc((*mol).r[i].na,sizeof(atom));
   l[i] = 0;
  }
  //Set the current residue, and residue index to 0
  rI = 0;
  curRes = ((*mol).r+0);
  for(i = 0; i < INWC; i++)
  {
    if(isAtom((ln+i)) == 1)
    {
//printf("in loop for isAtom\n");
      temp  = (*(ln+i)).f[8].c; sscanf(temp,"%d",&ntX);
      icode  = (*(ln+i)).f[9].c[0]; 
      chainid  = (*(ln+i)).f[7].c[0]; 
      if(((*curRes).n != ntX)||(icode!=curRes[0].IC[0])||(chainid!=curRes[0].cID[0]))//If the residue numbers do not match
      {//Seach through all of the resiudes until one is found that matches
       for(j = 0; j < (*mol).nr; j++){//then assign it to curRes
        if(mol[0].r[j].n == ntX){rI = j; break;}
       }
       curRes = ((*mol).r+rI);
      }
      curAtm = ((*curRes).a+l[rI]);
/*printf("\t still in loop for isAtom\n"); */
      //Getting the record name
      temp = (*(ln+i)).f[0].c; sscanf(temp,"%s",name);
      if(name!=NULL){(*curAtm).D = strdup(name);}
      else{(*curAtm).D = strdup("");}
      //Getting the atom #
      temp  = (*(ln+i)).f[1].c; sscanf(temp,"%d",&ntX);
      (*curAtm).n = ntX;
      //Getting the atom name
      temp  = (*(ln+i)).f[3].c; sscanf(temp,"%s",name);
      if(name!=NULL){(*curAtm).N = strdup(name);}
      else{(*curAtm).D = strdup("    ");}

      //Getting the chain identifier
      if((*(ln+i)).f[7].c!=NULL){
        (*curAtm).cID = (char*)calloc(2,sizeof(char));
        (*curAtm).cID[0] = (*(ln+i)).f[7].c[0];
        (*curAtm).cID[1] = '\0';
        }
      else{(*curAtm).cID = strdup(" ");}
/*printf("\t the (*curAtm).cID is >>>%s<<<\n",(*curAtm).cID); */

      //Getting the X coordinate
      temp  = (*(ln+i)).f[11].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.i = dblY;
      //Getting the Y coordinate
      temp  = (*(ln+i)).f[12].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.j = dblY;
      //Getting the Z coordinate
      temp  = (*(ln+i)).f[13].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.k = dblY;   

      l[rI]++;
//printf("the current atom's xyz are: %f  %f  %f \n",(*curAtm).x.i,(*curAtm).x.j,(*curAtm).x.k);
    }
  }

//printf("the number of atoms in the molecule's first residue is %d\n",(*mol).r[0].na);

 return mol;
}

