/*Author: Michael Tessier*/
#include <load_pdb.h>
//#include "../inc/load_pdb.h"
assembly* load_pdb(char* file_name)
{
  char* temp;
  WARN=1;SILENT=1;INWC=0;DATE=0;UNCUT=1;UNx=1;UNy=1;UNz=1;
  temp = suf; strcpy(temp,".pdb");
  temp = sufc; strcpy(temp,"_change.txt"); ACT='0';
  DEBUG=-1;LASTRES=-1;LASTOKX=0;LASTOKY=0;LASTOKZ=0;
  UNCTOL=0;CRYX=0;CRYY=0;CRYZ=0;LASTX=0;LASTY=0;LASTZ=0;
  //int resNum;
  int ma;
// int mi;
  //molecule* mol;
  assembly* asmbl;
  IN = myfopen(file_name, "r");
 /*scans pdb file for number of lines
  allocates memory for line structures
  initializes structures that contain line formats*/
  init_struct();
  rewind(IN);
  //Read in lines to ln structure
  printf("Reading in %s...\n",file_name);
  for(ma=0;ma<INWC;ma++){
  if(DEBUG>=1){printf("made it to here main line loop %d ...\n", ma);} 
  	rwm_line(ma+1); 
	}

  //Determine the number of molecules in the pdb
  printf("There are %d molecule(s)\n",howManyMolecules());
  asmbl = getAssembly();
  printf("PDB Information Successfully Read.\n");
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
 return (*getMolecule());
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
		counter ++;
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
  int i = 0;int j = 0;int k = 0;int l = 0;int atom_num;
  int init=0, next_mol_line=0;
  char in_molecule_switch = 'n';
  double x,y,z;
  char name[6] = ""; char* atmName = name;
  char atmElem[3];
  assembly* asmbl = (assembly*) calloc (1 , sizeof(assembly));
  (*asmbl).nm = molNum;
  molecule *curMol; residue* curRes;// residue* res;
  atom* curAtm;// atom* atm;
  /* made further changes for molecule double pointers in assembly
	starting on 20100723.  BLFoley */ 

  (*asmbl).m = (molecule**)calloc(molNum,sizeof(molecule*));//added by MNT on 20080806
  for(init=0;init<molNum;init++){ 
	(*asmbl).m[init] = (molecule*) calloc (1,sizeof(molecule)); 
	(*asmbl).m[init][0].nr=-1; /* initializing as empty/not-seen. */ 
	}
  i=0; 
  //Get all the information for the first molecule
/*  curMol = &(*asmbl).m[0][0];
  (*curMol).nr = findTotalResidue(i);
  (*curMol).r = (residue*) calloc ((*curMol).nr,sizeof(residue)); */
//printf("curMol nr is %d \n",(*curMol).nr);
  //Initialize the residue numbers to -1 so the program knows they're not used yet
/*  for(init = 0; init < (*curMol).nr; init++){(*curMol).r[init].n = -1;}
  i=0;
  next_mol_line=getResInfo((*curMol).r,i);
  curRes = ((*curMol).r+k);
  (*curRes).a = (atom*) calloc ((*curRes).na,sizeof(atom));
  curAtm = ((*curRes).a+l); */

//printf("First allocate of residues and atoms:\n");
//printf("\t molecule %d, residue %d -- nr is %d and na is %d\n",j,k,(*curMol).nr,(*curRes).na);

/* 
	i = line number
	j = molecule number
	k = residue number (within the molecule, not absolute)
	l = atom number ( I think... )
*/

in_molecule_switch = 'n';
for(i = 0; i < INWC; i++) {
//printf("entering the loop, is is %d -- the card is %s\n",i,(*(ln+i)).f[0].c); 
	if(j==molNum){break;}
	if(in_molecule_switch == 'n'){ /* If we are not in a molecule. */
		while( (i<(INWC)) && (isAtom(ln+i)==0) ){ i++; } /* advance to first ATOM/HETATM line */
		if(i==INWC){break;}
//printf("about to call getResInfo. i is %d -- j is %d\n",i,j);
		curMol = &asmbl[0].m[j][0];
		(*curMol).nr = findTotalResidue(i);
		(*curMol).r = (residue*) calloc ((*curMol).nr,sizeof(residue));
printf("\tcurMol nr is %d \n",(*curMol).nr);
		//Initialize the residue numbers to -1 so the program knows they're not used yet
		for(init = 0; init < (*curMol).nr; init++){(*curMol).r[init].n = -1;}
		next_mol_line=getResInfo((*curMol).r,i);
printf("\tcalled getResInfo. (*curMol).r[0].N is %s ; i is %d ; k is %d\n",(*curMol).r[0].N,i,k);
		k = 0;l = 0;
		curRes = &curMol[0].r[k];
printf("allocating residues and atoms:\n");
printf("\t molecule %d, residue %d -- nr is %d and na is %d\n",j,k,(*curMol).nr,(*curRes).na);
		(*curRes).a = (atom*) calloc ((*curRes).na,sizeof(atom));
		curAtm = ((*curRes).a+l);
		in_molecule_switch = 'y';/* At this point, we should be at an ATOM/HETATM entry. */
		j++;
		} 
	if(i==INWC){ /* If we ran out of lines... */
		if((*curMol).nr==-1){ /* If this is molecule does not currently contain residues */
			(*curMol).nr = 0;
			(*curMol).r = NULL;
			(*curMol).N = strdup("EMPTY_MOLECULE");
			if(j!=molNum+1){printf("Found %d molecules, but ran out of atoms before molecule %d.\n",molNum,j+1);
				printf("This is probably a very, very bad thing, but we're ignoring it for now.\n");}
			} 
		else if ((*curMol).nr>0){ /* if this molecule already contains residues */
			printf("Got to i=INWC while in a molecule.  Go fix code.\n"); }
		else { printf("Got to i=INWC for an uncoded value of nr.  Go fix code.\n");}
		}

	while( (i<INWC) && (in_molecule_switch=='y') ) 
	//while( (i<INWC) && (j<molNum) && (in_molecule_switch=='y') ) 
		{ /* while we are in a molecule */
		if(endOfMol((ln+i)) == 1){
			in_molecule_switch='n';
			if(next_mol_line!=i) {printf("i is %d and next_mol_line is %d ... should match.\n",i,next_mol_line);}
			break;
			}
		else {
			if(isAtom((ln+i)) == 1) {
				if(l == (*curRes).na) {
					k++;l = 0;
					curRes = ((*curMol).r+k);
//printf("allocating atoms only, molecule %d, residue %d and na is %d\n",j,k,(*curRes).na);
					(*curRes).a = (atom*) calloc ((*curRes).na,sizeof(atom));
					}
				curAtm = ((*curRes).a+l);
				//Getting the atom # //
				sscanf( (*(ln+i)).f[1].c ,"%d",&atom_num);
				//Getting the atom name //temp  = (*(ln+i)).f[3].c;
				sscanf( (*(ln+i)).f[3].c ,"%s",atmName);
				sscanf( (*(ln+i)).f[18].c ,"%s",atmElem); 
				//Getting the X, Y and Z coordinates
				sscanf( (*(ln+i)).f[11].c ,"%lf",&x);
				sscanf( (*(ln+i)).f[12].c ,"%lf",&y);
				sscanf( (*(ln+i)).f[13].c ,"%lf",&z);   
				(*curAtm).n = atom_num;//set the atom #
				(*curAtm).N = strdup(atmName);//set the atom name
				(*curAtm).E = strdup(atmElem);//set the atom element
				//set the atom's coordinates
				(*curAtm).x.i = x; (*curAtm).x.j = y; (*curAtm).x.k = z;
				l++;
				}
			}
		i++;
		} /* close while we are in a molecule */
	}
return asmbl;
}

int findTotalResidue(int start)
{
  int count=0;int i = start;int thisRes;//int curRes = 0;
  int inList,j;
  char* temp;
  int* known = (int*)calloc(count,sizeof(int));

//printf(" i is %d, and INWC is %d and ln[i].a-b are %d-%d\n",i,INWC,ln[i].a,ln[i].b);

  while( (i!=INWC) && (isAtom((ln+i))==0) ) { i++; }

  while( (i != INWC) && (endOfMol((ln+i)) == 0) )
  {
    if(isAtom((ln+i)) == 1)
    {
     temp  = (*(ln+i)).f[8].c;
     sscanf(temp,"%d",&thisRes);
     inList = 0;
     //The New Way
     for(j = 0; j < count; j++){
      if(known[j]==thisRes){inList=1; break;}
     }
     if(!inList){
      count++;
      known = (int*)realloc(known,count*sizeof(int));
      known[count-1] = thisRes;
     }
    }
   i++;
  }
  free(known);
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

//printf("**1.  i is %d; atom name is %s; atom number is %s\n",i,(*(ln+i)).f[3].c,(*(ln+i)).f[1].c);

while( (current_status == 0) && (i<INWC) ){ 
	current_status = isAtom(ln+i); /* find out if we have an atom line  */
	i++; } /* if not, check the next line */
if(i==INWC){ /* if we got to the end with no atoms */
	(*curRes).N = strdup("EMPTY_RESIDUE");
	(*curRes).n = -100;
	(*curRes).na = 0;
	return INWC;
	} 
i--; /* still here? decrement the counter to undo the while loop above */

//while(isAtom(ln+i) == 0 && i < INWC) 
//{
//i++;

while( (i != INWC)  && (endOfMol((ln+i)) == 0) ){ 
	if(isAtom(ln+i) == 1){
		temp  = (*(ln+i)).f[8].c; sscanf(temp,"%d",&resNum);
//printf("resNum is %d and (*curRes).n is %d\n",resNum,(*curRes).n);
//printf("**2.  i is %d; atom name is %s; atom number is %s\n",i,(*(ln+i)).f[3].c,(*(ln+i)).f[1].c);
		//If this is a different residue than the one in the previous line
		if(resNum != (*curRes).n){
			curRes = NULL;
			for(j = 0; j < count; j++){//Cycle through all of the found residues
				if((*(res+j)).n == resNum){ /* if already seen... */
					curRes = (res+j); break;} /* reset residue pointer */
				} 
			/* If it is not a known residue, assign it a new slot */
			if(curRes == NULL || curRes == 0x0){curRes = (res+count); count++;}
			temp = (*(ln+i)).f[5].c; sscanf(temp,"%s",resName);
			}
//printf("  -->  resNum is %d and (*curRes).n is %d\n",resNum,(*curRes).n);
//printf("**3.  i is %d; atom name is %s; atom number is %s\n",i,(*(ln+i)).f[3].c,(*(ln+i)).f[1].c);
		if((*curRes).n < 0){ /* If this residue has not been found yet */
			(*curRes).N = strdup(resName);	/* Set the residue name */
			(*curRes).n = resNum;		/* Set the actual residue number */
			(*curRes).na = 0;			/* ...and make sure the total # of atoms is 0 */
			} /* Otherwise we can assume all of these have already been set */
		(*curRes).na++;
//printf(" --> (curRes). n=%d, na=%d, N=%s\n",(*curRes).n,(*curRes).na,(*curRes).N);
//printf("**4.  i is %d; atom name is %s; atom number is %s\n",i,(*(ln+i)).f[3].c,(*(ln+i)).f[1].c);
		}//End if an atom 
	i++; //..and then incriment to the next line
	} //End loop through file
//} //end loop through file
return i;
}

molecule* getMolecule(void){
  int molNum = 1;
  int i = 0;//int j = 0;int k = 0;
  int ntX,j,rI,next_mol_line;
  double dblY;
  char* temp;char name[6] = ""; //char* atmName = name;
  molecule* mol = (molecule*)calloc(molNum,sizeof(molecule));
  residue* curRes;// residue* res;
  atom* curAtm;// atom* atm;
  //Get all the information for the first molecule
  (*mol).nr = findTotalResidue(i);
  (*mol).r = (residue*)calloc((*mol).nr,sizeof(residue));
  //Initialize the residue numbers to -1 so the program knows they're not used yet
  for(j = 0; j < (*mol).nr; j++){(*mol).r[j].n = -1;}
  next_mol_line=getResInfo((*mol).r,i);
  int l[(*mol).nr];
  for(i = 0; i < (*mol).nr; i++){
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
      temp  = (*(ln+i)).f[8].c; sscanf(temp,"%d",&ntX);
      if((*curRes).n != ntX)//If the residue numbers do not match
      {//Seach through all of the resiudes until one is found that matches
       for(j = 0; j < (*mol).nr; j++){//then assign it to curRes
        if(mol[0].r[j].n == ntX){rI = j; break;}
       }
       curRes = ((*mol).r+rI);
      }
      curAtm = ((*curRes).a+l[rI]);
      //Getting the record name
      temp = (*(ln+i)).f[0].c; sscanf(temp,"%s",name);
      (*curAtm).D = strdup(name);
      //Getting the atom #
      temp  = (*(ln+i)).f[1].c; sscanf(temp,"%d",&ntX);
      (*curAtm).n = ntX;
      //Getting the atom name
      temp  = (*(ln+i)).f[3].c; sscanf(temp,"%s",name);
      (*curAtm).N = strdup(name);
/*
These lines cause a seg fault for some reason.  They were removed in a 
previous commit, so I'm re-removing them...  20090819 BLF.
printf("Got atom number %d and name %s in get Molecule\n",(*curAtm).n,(*curAtm).N);
printf("in load_pdb_from_slurp: Atom number %d of residue number %d (%s) has the name %s\n",mol[0].r[5].a[21].n,mol[0].r[5].n,mol[0].r[5].N,mol[0].r[5].a[21].N);
*/
      //Getting the chain identifier
      (*curAtm).cID = (*(ln+i)).f[7].c[0];

      //Getting the X coordinate
      temp  = (*(ln+i)).f[11].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.i = dblY;
      //Getting the Y coordinate
      temp  = (*(ln+i)).f[12].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.j = dblY;
      //Getting the Z coordinate
      temp  = (*(ln+i)).f[13].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.k = dblY;   

      l[rI]++;
    }
  }

 return mol;
}

