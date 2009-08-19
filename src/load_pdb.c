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

 for(i = 0; i < INWC; i++)
  if(endOfMol((ln+i)) == 1)//If it is a TER, LINK, or CONECT card
   counter++;
 if(counter == 0 || isAtom(ln+(INWC-1)) == 1)
  counter ++;
 return counter;
}

assembly* getAssembly()
{
  int molNum = howManyMolecules();
  int i = 0;int j = 0;int k = 0;int l = 0;int atom_num;double x,y,z;
  char* temp;char name[6] = ""; char* atmName = name;
  assembly* asmbl = (assembly*) malloc (1 * sizeof(assembly));
  (*asmbl).nm = molNum;
  molecule* mol = (molecule*) malloc (molNum * sizeof(molecule));
  (*asmbl).m = (molecule**)calloc(1,sizeof(molecule*));//added by MNT on 20080806
  (*asmbl).m[0] = mol; // changed by BLF on 20080622 -- might need revisiting
  molecule* curMol; residue* curRes;// residue* res;
  atom* curAtm;// atom* atm;
  //Get all the information for the first molecule
  curMol = (mol+j);
  (*curMol).nr = findTotalResidue(i);
  (*curMol).r = (residue*) malloc ((*curMol).nr*sizeof(residue));
  getResInfo((*curMol).r,i);
  curRes = ((*curMol).r+k);
  (*curRes).a = (atom*) malloc ((*curRes).na*sizeof(atom));
  curAtm = ((*curRes).a+l);
  for(i = 0; i < INWC; i++)
  {
    //If it is a TER, LINK, or CONECT card
    if(endOfMol((ln+i)) == 1)
    {
      j++;
      curMol = (mol+j);
      (*curMol).nr = findTotalResidue(i+1);
      (*curMol).r = (residue*) malloc ((*curMol).nr*sizeof(residue));
      getResInfo((*curMol).r,i+1);
      curRes = ((*curMol).r+k);
      k = 0;l = 0;
      (*curRes).a = (atom*) malloc ((*curRes).na*sizeof(atom));
      curAtm = ((*curRes).a+l);
    }
    else if(isAtom((ln+i)) == 1)
    {
      if(l == (*curRes).na)
      {
       k++;l = 0;
       curRes = ((*curMol).r+k);
       (*curRes).a = (atom*) malloc ((*curRes).na*sizeof(atom)); 
       //curAtm = ((*curRes).a+l);
      }
      curAtm = ((*curRes).a+l);
      //Getting the atom #
      temp  = (*(ln+i)).f[1].c;
      sscanf(temp,"%d",&atom_num);
      //Getting the atom name
      temp  = (*(ln+i)).f[3].c;
      sscanf(temp,"%s",atmName);
      //Getting the X coordinate
      temp  = (*(ln+i)).f[11].c;
      sscanf(temp,"%lf",&x);
      //Getting the Y coordinate
      temp  = (*(ln+i)).f[12].c;
      sscanf(temp,"%lf",&y);
      //Getting the Z coordinate
      temp  = (*(ln+i)).f[13].c;
      sscanf(temp,"%lf",&z);   

      (*curAtm).n = atom_num;//set the atom #
      (*curAtm).N = strdup(atmName);//set the atom name
      //set the atom's coordinates
      (*curAtm).x.i = x; (*curAtm).x.j = y; (*curAtm).x.k = z;
      l++;
    }
  }
 return asmbl;
}

int findTotalResidue(int start)
{
  int count=0;int i = start;int thisRes;//int curRes = 0;
  int inList,j;
  char* temp;
  int* known = (int*)calloc(count,sizeof(int));
  while(endOfMol((ln+i)) == 0 && i != INWC)
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
  if(((*line).a == 4 && (*line).b == 2)||
     ((*line).a == 2 && ((*line).b == 9||(*line).b == 3)))
    return 1;
  else
    return 0;
}

int isAtom(linedef* line)
{
  if(((*line).a == 2 || (*line).a == 3) && (*line).b == 1)
    return 1;
  else
    return 0;
}
void getResInfo(residue* res, int start)
{
  int resNum,j;int i = start;
  int count = 0;
  residue* curRes = (res+0);
  char name [10];char* temp;char* resName = name;
  while(isAtom(ln+i) == 0 && i < INWC)
   i++;
  while(endOfMol((ln+i)) == 0 && i != INWC){
   if(isAtom(ln+i) == 1){
    temp  = (*(ln+i)).f[8].c; sscanf(temp,"%d",&resNum);
    //If this is a different residue than the one in the previous line
    if(resNum != (*curRes).n){
     curRes = NULL;
     for(j = 0; j < count; j++){//Cycle through all of the found residues
      if((*(res+j)).n == resNum){curRes = (res+j); break;}
     }				//If it is not a known residue, assign it a new slot
     if(curRes == NULL && curRes == 0x0){curRes = (res+count); count++;}
     temp = (*(ln+i)).f[5].c; sscanf(temp,"%s",resName);
    }
    //If this residue has not been found yet
    if((*curRes).n < 0){
     (*curRes).N = strdup(resName);	//Set the residue name
     (*curRes).n = resNum;		//Set the actual residue number
     (*curRes).na = 0;			//...and make sure the total # of atoms is 0
    }//Otherwise we can assume all of these have already been set
    (*curRes).na++;
   }//End if an atom
    //..and then incriment to the next line
    i++;
  }//End loop through file
 return ;
}

molecule* getMolecule(void){
  int molNum = 1;
  int i = 0;//int j = 0;int k = 0;
  int ntX,j,rI;
  double dblY;
  char* temp;char name[6] = ""; //char* atmName = name;
  molecule* mol = (molecule*)calloc(molNum,sizeof(molecule));
  residue* curRes;// residue* res;
  atom* curAtm;// atom* atm;
  //Get all the information for the first molecule
  (*mol).nr = findTotalResidue(i);
  (*mol).r = (residue*)calloc((*mol).nr,sizeof(residue));
  //Initialize of the residue numbers to -1 so the program knows they're not used yet
  for(j = 0; j < (*mol).nr; j++){(*mol).r[j].n = -1;}
  getResInfo((*mol).r,i);
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

